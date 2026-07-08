// ESM-C (EvolutionaryScale) sequence-pLM tier.
// Embed → kNN + nearest-centroid retrieval → contrastive head.
// Real execution needs a GPU + the project venv on met
// (/array1/xinpeng/dbcan4-advanced/venv). The stub block reproduces the exact
// published layout with no GPU so the DAG and the downstream contract can be
// proven anywhere (-stub-run).

process ESMC_EMBED {
    tag   "${sample}"
    label 'gpu'
    conda "${params.esmc_env}"

    input:
    tuple val(sample), path(faa)

    output:
    tuple val(sample), path("${sample}.esmc.npz"), emit: emb

    script:
    """
    embed_esmc.py --fasta ${faa} --out-prefix ${sample}.esmc \\
        --model ${params.esmc_model} --nshards 1 --shard 0
    # single-shard runs write <prefix>.shard0.npz; expose the canonical name
    mv ${sample}.esmc.shard0.npz ${sample}.esmc.npz
    """

    stub:
    """
    : > ${sample}.esmc.npz
    """
}

// kNN + nearest-centroid share one retrieval script (two schemes, one pass).
// It writes a wide per-query TSV; the pipeline then normalizes each scheme's
// columns into a standard §2.1 file (ESM-C-kNN.tsv, ESM-C-centroid.tsv).
process ESMC_RETRIEVAL {
    tag   "${sample}"
    label 'cpu'
    conda "${params.esmc_env}"
    publishDir { "${params.outdir}/cazyme_advanced/predictions/${sample}" }, mode: 'copy', pattern: '*.tsv'
    input:
    tuple val(sample), path(emb)
    path ref_emb_dir
    path eval_labels

    output:
    tuple val(sample), path("ESM-C-kNN.tsv"), path("ESM-C-centroid.tsv"), emit: preds
    path "esmc_retrieval_raw.tsv", emit: raw

    script:
    """
    retrieval_esmc.py \\
        --ref-prefix ${ref_emb_dir}/reference_2024.esmc \\
        --eval-prefix ${sample}.esmc \\
        --labels ${eval_labels} \\
        --out-summary esmc_retrieval_summary.json \\
        --out-pred esmc_retrieval_raw.tsv --k ${params.knn_k}
    normalize_predictions.py --tool ESM-C-kNN      --in esmc_retrieval_raw.tsv --out ESM-C-kNN.tsv
    normalize_predictions.py --tool ESM-C-centroid --in esmc_retrieval_raw.tsv --out ESM-C-centroid.tsv
    """

    stub:
    """
    cp ${projectDir}/assets/stub/esmc_retrieval_raw.tsv esmc_retrieval_raw.tsv
    normalize_predictions.py --tool ESM-C-kNN      --in esmc_retrieval_raw.tsv --out ESM-C-kNN.tsv
    normalize_predictions.py --tool ESM-C-centroid --in esmc_retrieval_raw.tsv --out ESM-C-centroid.tsv
    """
}

process ESMC_CONTRASTIVE {
    tag   "${sample}"
    label 'gpu'
    conda "${params.esmc_env}"
    publishDir { "${params.outdir}/cazyme_advanced/predictions/${sample}" }, mode: 'copy', pattern: 'ESM-C-contrastive.tsv'
    input:
    tuple val(sample), path(emb)
    path ref_emb_dir
    path eval_labels

    output:
    tuple val(sample), path("ESM-C-contrastive.tsv"), emit: preds

    script:
    """
    train_heads.py \\
        --ref-prefix ${ref_emb_dir}/reference_2024.esmc \\
        --eval-prefix ${sample}.esmc \\
        --labels ${eval_labels} \\
        --out-pred head_raw.tsv --out-metrics head_metrics.json
    normalize_predictions.py --tool ESM-C-contrastive --in head_raw.tsv --out ESM-C-contrastive.tsv
    """

    stub:
    """
    cp ${projectDir}/assets/stub/head_raw.tsv head_raw.tsv
    normalize_predictions.py --tool ESM-C-contrastive --in head_raw.tsv --out ESM-C-contrastive.tsv
    """
}
