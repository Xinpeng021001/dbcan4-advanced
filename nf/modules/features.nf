// Per-protein feature tier (deeper analysis of a candidate CAZyme):
// SignalP 6.0 (secretion), DeepTMHMM (membrane topology). Both are the standard
// tools; here they publish the contract's §2.2/§2.3 TSVs. Stubs ship
// contract-shaped rows so the layout is proven without the tools installed.

process SIGNALP6 {
    tag   "${sample}"
    label 'gpu'
    conda "${params.signalp_env}"
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'signalp6.tsv'
    input:
    tuple val(sample), path(faa)

    output:
    tuple val(sample), path("signalp6.tsv"), emit: feat

    script:
    """
    signalp6 --fastafile ${faa} --organism eukarya --format none \\
        --output_dir signalp_out --mode fast
    signalp6_to_tsv.py --in signalp_out/prediction_results.txt --out signalp6.tsv
    """

    stub:
    """
    cp ${projectDir}/assets/stub/signalp6.tsv signalp6.tsv
    """
}

process DEEPTMHMM {
    tag   "${sample}"
    label 'gpu'
    conda "${params.deeptmhmm_env}"
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'deeptmhmm.tsv'
    input:
    tuple val(sample), path(faa)

    output:
    tuple val(sample), path("deeptmhmm.tsv"), emit: feat

    script:
    """
    biolib run DTU/DeepTMHMM --fasta ${faa}
    deeptmhmm_to_tsv.py --in biolib_results/TMRs.gff3 --out deeptmhmm.tsv
    """

    stub:
    """
    cp ${projectDir}/assets/stub/deeptmhmm.tsv deeptmhmm.tsv
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
    write_manifest.py --outdir ${params.outdir} --sample ${sample} \\
        --release-label '${params.release_label}' \\
        --pipeline-version '${params.pipeline_version}' \\
        --tool-versions '${params.tool_versions_json}'
    cp ${params.outdir}/cazyme_advanced/manifest.json manifest.json
    """
}
