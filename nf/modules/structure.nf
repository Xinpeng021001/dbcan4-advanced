// Structure tier — gated onto hard/baseline-missed proteins.
// (a) ESMFold the subset → real 3Di.
// (b) Foldseek 3Di search vs the CAZyme3D structure DB (family from best hit).
// (c) SaProt structure-aware embedding retrieval (orthogonal structure signal).
// Real execution needs GPU (ESMFold) + Foldseek + the CAZyme3D DB on met. The
// stubs publish contract-shaped outputs (incl. a served PDB) with no compute.

process ESMFOLD {
    tag   "${sample}"
    label 'gpu'
    conda "${params.esmfold_env}"
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'structures*'
    input:
    tuple val(sample), path(faa)

    output:
    tuple val(sample), path("structures.tsv"), path("structures/"), emit: struct

    script:
    """
    mkdir -p structures
    esmfold_infer.py --fasta ${faa} --out-dir structures --index structures.tsv \\
        --max-proteins ${params.struct_max}
    """

    stub:
    """
    mkdir -p structures
    cp ${projectDir}/assets/stub/structures/*.pdb structures/ 2>/dev/null || true
    cp ${projectDir}/assets/stub/structures.tsv structures.tsv
    """
}

process FOLDSEEK_CAZYME3D {
    tag   "${sample}"
    label 'cpu'
    conda "${params.foldseek_env}"
    publishDir { "${params.outdir}/cazyme_advanced/predictions/${sample}" }, mode: 'copy', pattern: 'Foldseek-CAZyme3D.tsv'
    input:
    tuple val(sample), path(struct_tsv), path(struct_dir)

    output:
    tuple val(sample), path("Foldseek-CAZyme3D.tsv"), emit: preds

    script:
    """
    foldseek easy-search ${struct_dir} ${params.cazyme3d_db} foldseek_raw.m8 tmp \\
        --format-output 'query,target,fident,alntmscore,lddt,prob' -e ${params.foldseek_evalue}
    foldseek_to_family.py --hits foldseek_raw.m8 --db-labels ${params.cazyme3d_labels} \\
        --out foldseek_raw.tsv
    normalize_predictions.py --tool Foldseek-CAZyme3D --in foldseek_raw.tsv --out Foldseek-CAZyme3D.tsv
    """

    stub:
    """
    cp ${projectDir}/assets/stub/foldseek_raw.tsv foldseek_raw.tsv
    normalize_predictions.py --tool Foldseek-CAZyme3D --in foldseek_raw.tsv --out Foldseek-CAZyme3D.tsv
    """
}

process SAPROT {
    tag   "${sample}"
    label 'gpu'
    conda "${params.saprot_env}"
    publishDir { "${params.outdir}/cazyme_advanced/predictions/${sample}" }, mode: 'copy', pattern: 'SaProt.tsv'
    input:
    tuple val(sample), path(struct_tsv), path(struct_dir)
    path ref_saprot_dir

    output:
    tuple val(sample), path("SaProt.tsv"), emit: preds

    script:
    """
    saprot_retrieval.py --struct-index ${struct_tsv} --struct-dir ${struct_dir} \\
        --ref-prefix ${ref_saprot_dir}/reference_2024.saprot --out saprot_raw.tsv
    normalize_predictions.py --tool SaProt --in saprot_raw.tsv --out SaProt.tsv
    """

    stub:
    """
    cp ${projectDir}/assets/stub/saprot_raw.tsv saprot_raw.tsv
    normalize_predictions.py --tool SaProt --in saprot_raw.tsv --out SaProt.tsv
    """
}
