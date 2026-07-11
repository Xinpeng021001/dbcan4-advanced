// Per-protein feature tier (deeper analysis of a candidate CAZyme):
// membrane topology + secretion (DeepTMHMM, real BioLib-cloud run) and signal
// peptide. DeepTMHMM is the real workhorse here: it predicts BOTH the
// transmembrane topology (§2.3 tm_topology) AND the N-terminal signal peptide,
// so a single real run populates deeptmhmm.tsv and, honestly, the signal_peptide
// TSV when the licensed SignalP-6.0 is not installed on the host.
//
// Each process has a real script: block (venv absolute paths — conda directives
// are not reliable on met) and a stub: block (canned contract-shaped rows) so the
// DAG proves anywhere without the tools.

// DeepTMHMM — real membrane topology + signal peptide via pybiolib -> BioLib cloud.
// GOTCHA: `biolib run --fasta X` stages the arg file by BASENAME on the cloud, so
// we copy the FASTA to the task dir and pass a bare filename (an absolute path
// fails cloud-side with FileNotFoundError). Emits BOTH deeptmhmm.tsv (tm_topology)
// and signalp6_from_dtm.tsv (signal_peptide, derived) from one run.
process DEEPTMHMM {
    tag   "${sample}"
    label 'cpu'
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'deeptmhmm.tsv'
    input:
    tuple val(sample), path(faa)

    output:
    tuple val(sample), path("deeptmhmm.tsv"),         emit: feat
    tuple val(sample), path("signalp6_from_dtm.tsv"), emit: sp_derived

    script:
    """
    cp ${faa} dtm_input.faa
    ${params.biolib_bin} run DTU/DeepTMHMM --fasta dtm_input.faa
    deeptmhmm_to_tsv.py --gff3 biolib_results/TMRs.gff3 \\
        --three-line biolib_results/predicted_topologies.3line \\
        --out-tm deeptmhmm.tsv --out-sp signalp6_from_dtm.tsv
    """

    stub:
    """
    cp ${projectDir}/assets/stub/deeptmhmm.tsv deeptmhmm.tsv
    cp ${projectDir}/assets/stub/signalp6.tsv  signalp6_from_dtm.tsv
    """
}

// SignalP-6.0 (secretion) — licensed, may not be installed. HONEST FALLBACK:
// if the signalp6 binary is on PATH we run it; otherwise we publish the
// DeepTMHMM-derived signal_peptide TSV (clearly sourced in its `extra` column).
// Never fabricates an SP probability.
process SIGNALP6 {
    tag   "${sample}"
    label 'cpu'
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'signalp6.tsv'
    input:
    tuple val(sample), path(faa), path(sp_derived)

    output:
    tuple val(sample), path("signalp6.tsv"), emit: feat

    script:
    """
    if command -v signalp6 >/dev/null 2>&1; then
        signalp6 --fastafile ${faa} --organism eukarya --format none \\
            --output_dir signalp_out --mode fast
        signalp6_to_tsv.py --in signalp_out/prediction_results.txt --out signalp6.tsv
    else
        echo "[SIGNALP6] signalp6 not on PATH — using DeepTMHMM-derived signal peptide (honest fallback)" >&2
        cp ${sp_derived} signalp6.tsv
    fi
    """

    stub:
    """
    cp ${projectDir}/assets/stub/signalp6.tsv signalp6.tsv
    """
}

// Fusion — combine sequence-pLM + structure evidence into one calibrated call.
// Agreement across orthogonal signals is the strongest remote-homolog evidence.
process FUSION {
    tag   "${sample}"
    label 'cpu'
    conda "${params.fusion_env}"
    publishDir { "${params.outdir}/cazyme_advanced/predictions/${sample}" }, mode: 'copy', pattern: 'fusion.tsv'
    input:
    tuple val(sample), path(pred_tsvs)

    output:
    tuple val(sample), path("fusion.tsv"), emit: preds

    script:
    """
    fuse_predictions.py --inputs ${pred_tsvs} --out fusion_raw.tsv \\
        --weights ${params.fusion_weights} --min-confidence ${params.fusion_min_conf}
    normalize_predictions.py --tool fusion --in fusion_raw.tsv --out fusion.tsv
    """

    stub:
    """
    fuse_predictions.py --inputs ${pred_tsvs} --out fusion_raw.tsv \\
        --weights ${params.fusion_weights} --min-confidence ${params.fusion_min_conf}
    normalize_predictions.py --tool fusion --in fusion_raw.tsv --out fusion.tsv
    """
}

process COLLATE_MANIFEST {
    tag   "${sample}"
    label 'cpu'
    publishDir "${params.outdir}/cazyme_advanced", mode: 'copy', pattern: 'manifest.json'

    input:
    tuple val(sample), path(pred_files), path(feat_files)

    output:
    path "manifest.json", emit: manifest

    script:
    """
    write_manifest.py --outdir ${params.outdir} --sample ${sample} --stage-dir . \\
        --release-label '${params.release_label}' \\
        --pipeline-version '${params.pipeline_version}' \\
        --tool-versions '${params.tool_versions_json}'
    """
}
