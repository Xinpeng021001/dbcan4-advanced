"""Advanced ingester: loads the standard contract as a new release, attaching
calls to existing baseline genes; idempotent; flags advanced-only families."""
from __future__ import annotations

import json
from pathlib import Path

from sqlalchemy import func, select

from bioforge.ingest.loader import ingest_directory
from bioforge.ingest.loader_advanced import ingest_advanced_manifest
from bioforge.models import CazymeAnnotation, Gene, ProteinFeature, Release, Sample


def _write_advanced_fixture(tmp: Path) -> Path:
    """A minimal cazyme_advanced/ tree for sampleA (genes PROKKA_00001/2/3/5
    from the shared sample_data fixture). demo covers: agreement (00001),
    an ADVANCED-ONLY call on a gene the baseline also called, and features."""
    base = tmp / "cazyme_advanced"
    (base / "predictions" / "sampleA").mkdir(parents=True)
    (base / "features" / "sampleA").mkdir(parents=True)
    (base / "features" / "sampleA" / "structures").mkdir(parents=True)

    def tsv(path, header, rows):
        with open(path, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    hdr = ["protein_id", "family", "confidence", "ec", "all_families", "extra"]
    # ESM-C-kNN: recovers GH13 on 00001 (agrees w/ baseline) + AA9 on 00005
    tsv(base / "predictions/sampleA/ESM-C-kNN.tsv", hdr, [
        ["PROKKA_00001", "GH13", "0.98", "-", "-", "{}"],
        ["PROKKA_00005", "AA9", "0.91", "-", "-", '{"knn_purity":1.0}'],
    ])
    # fusion: consensus GH13 (00001) + AA9 (00005), the advanced-only family
    tsv(base / "predictions/sampleA/fusion.tsv", hdr, [
        ["PROKKA_00001", "GH13", "0.95", "-", "-", '{"agreement":2,"signals":["sequence"]}'],
        ["PROKKA_00005", "AA9", "0.88", "-", "-", '{"agreement":2,"signals":["sequence","structure"]}'],
    ])
    # SignalP6 feature
    tsv(base / "features/sampleA/signalp6.tsv",
        ["protein_id", "prediction", "sp_prob", "cs_position", "extra"],
        [["PROKKA_00001", "SP", "0.97", "20", "{}"],
         ["PROKKA_00002", "NO_SP", "0.02", "-", "{}"]])
    # DeepTMHMM feature
    tsv(base / "features/sampleA/deeptmhmm.tsv",
        ["protein_id", "prediction", "n_tm", "topology", "extra"],
        [["PROKKA_00001", "Globular", "0", "-", "{}"]])
    # structure feature + a PDB file
    (base / "features/sampleA/structures/PROKKA_00001.pdb").write_text(
        "REMARK test\nATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\nEND\n")
    tsv(base / "features/sampleA/structures.tsv",
        ["protein_id", "source", "plddt", "length", "path", "extra"],
        [["PROKKA_00001", "ESMFold", "88.5", "129", "structures/PROKKA_00001.pdb", "{}"]])

    manifest = {
        "contract_version": "1.0", "pipeline": "dbcan4-advanced",
        "pipeline_version": "0.1.0", "release_label": "advanced-test",
        "release_notes": "unit-test advanced calls",
        "tool_versions": {"esm": "3.2.1"},
        "samples": [{
            "sample_key": "sampleA",
            "cazyme_predictions": [
                {"tool": "ESM-C-kNN", "path": "predictions/sampleA/ESM-C-kNN.tsv"},
                {"tool": "fusion", "path": "predictions/sampleA/fusion.tsv"},
            ],
            "protein_features": [
                {"feature_type": "signal_peptide", "tool": "SignalP6",
                 "path": "features/sampleA/signalp6.tsv"},
                {"feature_type": "tm_topology", "tool": "DeepTMHMM",
                 "path": "features/sampleA/deeptmhmm.tsv"},
                {"feature_type": "structure", "tool": "ESMFold",
                 "path": "features/sampleA/structures.tsv"},
            ],
        }],
    }
    (base / "manifest.json").write_text(json.dumps(manifest))
    return base / "manifest.json"


def test_advanced_ingest_creates_new_release_and_attaches_calls(
        Session, sample_data, tmp_path):
    # Baseline first.
    with Session() as s:
        ingest_directory(s, sample_data, label="baseline")
    man = _write_advanced_fixture(tmp_path)

    with Session() as s:
        rep = ingest_advanced_manifest(s, man, structures_dir=tmp_path / "served")

    assert rep.release_id is not None
    assert rep.samples_matched == 1
    assert rep.advanced_calls_added == 4          # 2 kNN + 2 fusion
    assert rep.per_tool == {"ESM-C-kNN": 2, "fusion": 2}
    assert rep.features_added == 4                # 2 SP + 1 TM + 1 structure
    assert rep.structures_copied == 1
    assert not rep.samples_missing

    with Session() as s:
        # A NEW release exists; baseline sample was NOT duplicated.
        assert s.scalar(select(func.count(Release.id))) == 2
        assert s.scalar(select(func.count(Sample.id))) == 2  # sampleA + sampleB only
        # Advanced calls carry confidence + method_family='advanced'.
        adv = s.scalars(select(CazymeAnnotation).where(
            CazymeAnnotation.method_family == "advanced")).all()
        assert len(adv) == 4
        assert all(c.confidence is not None for c in adv)
        assert all(c.release_id == rep.release_id for c in adv)
        # AA9 is ADVANCED-ONLY: baseline never called it.
        base_fams = {c.cazy_family for c in s.scalars(select(CazymeAnnotation).where(
            CazymeAnnotation.method_family == "baseline")).all()}
        adv_fams = {c.cazy_family for c in adv}
        assert "AA9" in adv_fams and "AA9" not in base_fams
        # Features landed with a structure path.
        feats = s.scalars(select(ProteinFeature)).all()
        assert len(feats) == 4
        struct = [f for f in feats if f.feature_type == "structure"][0]
        assert struct.structure_path and struct.structure_path.endswith(".pdb")
        assert struct.score == 88.5


def test_advanced_ingest_is_idempotent(Session, sample_data, tmp_path):
    with Session() as s:
        ingest_directory(s, sample_data, label="baseline")
    man = _write_advanced_fixture(tmp_path)
    with Session() as s:
        ingest_advanced_manifest(s, man, structures_dir=tmp_path / "served")
    with Session() as s:
        before = s.scalar(select(func.count(CazymeAnnotation.id)))
        rel_before = s.scalar(select(func.count(Release.id)))
    with Session() as s:
        rep2 = ingest_advanced_manifest(s, man, structures_dir=tmp_path / "served")
    with Session() as s:
        after = s.scalar(select(func.count(CazymeAnnotation.id)))
        rel_after = s.scalar(select(func.count(Release.id)))
    assert rep2.release_id is None            # no-op
    assert before == after                    # no duplicate calls
    assert rel_before == rel_after            # no phantom release


def test_advanced_ingest_raw_column_mapping(Session, sample_data, tmp_path):
    """Raw/legacy TSV mode: map arbitrary columns via id_col/family_col."""
    with Session() as s:
        ingest_directory(s, sample_data, label="baseline")
    base = tmp_path / "cazyme_advanced"
    (base / "predictions" / "sampleA").mkdir(parents=True)
    # A prototype-style wide TSV (query_id/cent_pred/cent_conf).
    with open(base / "predictions/sampleA/raw.tsv", "w") as fh:
        fh.write("query_id\tnovelty\tcent_pred\tcent_conf\n")
        fh.write("PROKKA_00002\tnovel_seq\tGH2\t0.87\n")
    (base / "manifest.json").write_text(json.dumps({
        "contract_version": "1.0", "pipeline": "dbcan4-advanced",
        "release_label": "advanced-raw",
        "samples": [{"sample_key": "sampleA", "cazyme_predictions": [
            {"tool": "ESM-C-centroid", "path": "predictions/sampleA/raw.tsv",
             "id_col": "query_id", "family_col": "cent_pred",
             "confidence_col": "cent_conf"}]}],
    }))
    with Session() as s:
        rep = ingest_advanced_manifest(s, base / "manifest.json",
                                       structures_dir=tmp_path / "served")
    assert rep.advanced_calls_added == 1
    with Session() as s:
        c = s.scalars(select(CazymeAnnotation).where(
            CazymeAnnotation.tool == "ESM-C-centroid")).first()
        assert c is not None and c.cazy_family == "GH2"
        assert abs(c.confidence - 0.87) < 1e-6


def test_v11_feature_parsers(tmp_path):
    """OUTPUT_CONTRACT v1.1 feature types (domain, structure_hit, localization,
    physicochem, ec_prediction) all parse onto the generic ProteinFeature shape."""
    from bioforge.ingest.parse_advanced import read_features, FeatFile

    d = tmp_path
    (d / "domains.tsv").write_text(
        "protein_id\tacc\tname\tstart\tend\tevalue\tscore\textra\n"
        "P1\tPF00295\tGlyco_hydro_28\t62\t412\t3.7e-36\t125.6\t{\"hmm_coverage\":0.88}\n")
    (d / "hits.tsv").write_text(
        "protein_id\ttarget\ttarget_family\ttmscore\tprob\tlddt\tevalue\textra\n"
        "P1\tQKX54669.1\tGH28\t0.82\t1.0\t0.64\t2.9e-14\t{\"bits\":582}\n")
    (d / "ec.tsv").write_text(
        "protein_id\tec_number\tconfidence\trank\ttool\textra\n"
        "P1\t3.2.1.40\t0.1053\t1\tCLEAN\t{\"confidence_type\":\"gmm\"}\n")
    (d / "pc.tsv").write_text(
        "protein_id\tmw_da\tpi\tinstability\tgravy\taromaticity\textra\n"
        "P1\t117294\t5.18\t29.42\t-0.19\t0.106\t{}\n")
    (d / "loc.tsv").write_text(
        "protein_id\tlocalization\tconfidence\tmethod\textra\n"
        "P1\tExtracellular\t-\tderived\t{}\n")

    dom = list(read_features(FeatFile("domain", "Pfam/hmmscan", d / "domains.tsv")))
    assert dom[0].label == "Glyco_hydro_28" and dom[0].start == 62 and dom[0].end == 412
    assert dom[0].attributes["acc"] == "PF00295"

    hit = list(read_features(FeatFile("structure_hit", "Foldseek-CAZyme3D", d / "hits.tsv")))
    assert hit[0].label == "GH28" and abs(hit[0].score - 0.82) < 1e-6
    assert hit[0].attributes["target"] == "QKX54669.1"

    ec = list(read_features(FeatFile("ec_prediction", "CLEAN", d / "ec.tsv")))
    assert ec[0].label == "3.2.1.40" and abs(ec[0].score - 0.1053) < 1e-6

    pc = list(read_features(FeatFile("physicochem", "Biopython", d / "pc.tsv")))
    assert abs(pc[0].score - 117294) < 1 and abs(pc[0].attributes["pi"] - 5.18) < 1e-6

    loc = list(read_features(FeatFile("localization", "DeepLoc", d / "loc.tsv")))
    assert loc[0].label == "Extracellular" and loc[0].attributes["method"] == "derived"


def test_feature_tools_registry():
    """Every v1.1 feature tool is registered with a colour + status."""
    from bioforge import methods
    for key in ("Pfam/hmmscan", "ESMFold", "Foldseek-CAZyme3D", "DeepTMHMM",
                "SignalP6", "CLEAN", "Biopython", "DeepLoc"):
        ft = methods.feature_tool(key)
        assert ft is not None, key
        assert ft.colour.startswith("#")
        assert ft.status in ("real_tool", "scaffold")
    assert methods.feature_tool("SignalP6").status == "scaffold"
    assert methods.feature_tool("Pfam/hmmscan").status == "real_tool"
