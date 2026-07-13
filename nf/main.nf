#!/usr/bin/env nextflow
/*
 * dbCAN4-advanced — advanced fungal-protein CAZyme annotation pipeline
 * ====================================================================
 * Wraps the advanced tiers (ESM-C sequence-pLM, ESMFold+Foldseek/CAZyme3D,
 * SaProt, SignalP6, DeepTMHMM, fusion) as a reproducible Nextflow DSL2 workflow
 * that publishes the STANDARD OUTPUT CONTRACT (nf/OUTPUT_CONTRACT.md). BioForge
 * then ingests the published manifest as a new versioned release.
 *
 *   proteins.faa ─┬─ ESM-C embed ─┬─ kNN + centroid ─┐
 *                 │               └─ contrastive     │
 *                 ├─ ESMFold ─┬─ Foldseek/CAZyme3D ──┼─ FUSION ─┐
 *                 │           └─ SaProt ─────────────┘          │
 *                 ├─ SignalP6 ──────────────(feature)           ├─ manifest.json
 *                 └─ DeepTMHMM ─────────────(feature)           │
 *                                                     structures┘
 *
 * Run modes:
 *   -profile stub    prove the DAG + contract with no GPU/tools (CI, laptops)
 *   -profile met     real execution on met (GPU, conda envs, CAZyme3D DB)
 */
nextflow.enable.dsl = 2

include { ESMC_EMBED; ESMC_RETRIEVAL; ESMC_CONTRASTIVE } from './modules/esmc.nf'
include { ESMFOLD; FOLDSEEK_CAZYME3D; SAPROT }           from './modules/structure.nf'
include { SIGNALP6; DEEPTMHMM; FUSION; COLLATE_MANIFEST } from './modules/features.nf'
include { BASELINE_DBCAN }                                  from './modules/baseline.nf'
include { INTERPROSCAN }                                     from './modules/interproscan.nf'
include { PFAM_DOMAINS; STRUCTURE_HITS; LOCALIZATION; PHYSICOCHEM; CLEAN_EC } from './modules/features_extra.nf'

def helpMessage() {
    log.info """
    dbCAN4-advanced annotation pipeline
    Usage:
      nextflow run nf/main.nf -profile stub  --input samplesheet.csv --outdir results
      nextflow run nf/main.nf -profile met   --input samplesheet.csv --outdir results

    samplesheet.csv columns: sample,faa   (protein FASTA per sample)
    Key params (see nextflow.config): --esmc_env --foldseek_env --cazyme3d_db
                                      --struct_max --knn_k --release_label
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); return }

    // --- input: samplesheet (sample,faa) ---
    Channel.fromPath(params.input, checkIfExists: true)
        | splitCsv(header: true)
        | map { row -> tuple(row.sample, file(row.faa, checkIfExists: true)) }
        | set { ch_faa }

    // --- Baseline dbCAN (run_dbcan) -> funcscan layout for BioForge ---
    BASELINE_DBCAN(ch_faa)

    ref_emb   = file(params.ref_emb_dir)
    ref_sap   = file(params.ref_saprot_dir)
    ev_labels = file(params.eval_labels)

    // --- Tier 1: ESM-C sequence-pLM ---
    ESMC_EMBED(ch_faa)
    ESMC_RETRIEVAL(ESMC_EMBED.out.emb, ref_emb, ev_labels)
    ESMC_CONTRASTIVE(ESMC_EMBED.out.emb, ref_emb, ev_labels)

    // --- Tier 2: structure (ESMFold → Foldseek/CAZyme3D + SaProt) ---
    ESMFOLD(ch_faa)
    FOLDSEEK_CAZYME3D(ESMFOLD.out.struct)
    SAPROT(ESMFOLD.out.struct, ref_sap)

    // --- Per-protein features ---
    // DeepTMHMM runs first: it predicts BOTH tm topology and the signal peptide,
    // and its derived signal_peptide feeds SIGNALP6's honest fallback when the
    // licensed signalp6 binary is not installed.
    DEEPTMHMM(ch_faa)
    SIGNALP6(ch_faa.join(DEEPTMHMM.out.sp_derived))

    // --- Comprehensive v1.1 features (§2.5-2.9) ---
    PFAM_DOMAINS(ch_faa)
    // InterProScan (or the offline Pfam->GO fallback) -> funcscan tree.
    // Reuses the Pfam domains so GO + InterPro populate the web gene page.
    INTERPROSCAN(ch_faa.join(PFAM_DOMAINS.out.feat))
    PHYSICOCHEM(ch_faa)
    CLEAN_EC(ch_faa)
    LOCALIZATION(SIGNALP6.out.feat)
    STRUCTURE_HITS(ESMFOLD.out.struct)

    // --- Fusion: gather every per-method prediction TSV per sample ---
    ESMC_RETRIEVAL.out.preds
        .map { s, knn, cent -> tuple(s, [knn, cent]) }
        .join(ESMC_CONTRASTIVE.out.preds.map { s, f -> tuple(s, [f]) })
        .join(FOLDSEEK_CAZYME3D.out.preds.map { s, f -> tuple(s, [f]) })
        .join(SAPROT.out.preds.map { s, f -> tuple(s, [f]) })
        .map { s, a, b, c, d -> tuple(s, (a + b + c + d)) }
        .set { ch_all_preds }

    FUSION(ch_all_preds)

    // --- Collate the manifest once all published files exist for the sample ---
    // Bundle every prediction file (incl. fusion) + every feature file per sample.
    ch_all_preds
        .join(FUSION.out.preds.map { s, f -> tuple(s, [f]) })
        .map { s, preds, fus -> tuple(s, preds + fus) }
        .join(
            SIGNALP6.out.feat.map { s, f -> tuple(s, [f]) }
                .join(DEEPTMHMM.out.feat.map { s, f -> tuple(s, [f]) })
                .join(ESMFOLD.out.struct.map { s, tsv, dir -> tuple(s, [tsv]) })
                .join(PFAM_DOMAINS.out.feat.map { s, f -> tuple(s, [f]) })
                .join(STRUCTURE_HITS.out.feat.map { s, f -> tuple(s, [f]) })
                .join(LOCALIZATION.out.feat.map { s, f -> tuple(s, [f]) })
                .join(PHYSICOCHEM.out.feat.map { s, f -> tuple(s, [f]) })
                .join(CLEAN_EC.out.feat.map { s, f -> tuple(s, [f]) })
                .map { s, a, b, c, d, e, f, g, h -> tuple(s, a + b + c + d + e + f + g + h) }
        )
        .map { s, preds, feats -> tuple(s, preds, feats) }
        .set { ch_collate }

    COLLATE_MANIFEST(ch_collate)

    workflow.onComplete = {
        log.info (workflow.success
            ? "OK dbCAN4-advanced complete -> ${params.outdir}/cazyme_advanced/manifest.json"
            : "FAILED dbCAN4-advanced pipeline")
    }
}
