// Comprehensive per-protein feature tier (v1.1 §2.5-2.9) — the five features that
// power the hero deep-dive page: Pfam/hmmscan domains, Foldseek structure_hits,
// derived localization, Biopython physicochem, CLEAN EC prediction. Each process
// has a real script: block (the tool + a feature_converters.py subcommand) and a
// stub: block (canned contract-shaped asset) so the DAG proves anywhere.

process PFAM_DOMAINS {
    tag   "${sample}"
    label 'cpu'
    conda "${params.dbcan_env}"
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'domains.tsv'
    input:
    tuple val(sample), path(faa)
    output:
    tuple val(sample), path("domains.tsv"), emit: feat
    script:
    """
    hmmscan --domtblout domains.domtbl --cut_ga -o /dev/null ${params.pfam_hmm} ${faa}
    feature_converters.py domains --domtbl domains.domtbl --out domains.tsv
    """
    stub:
    """
    cp ${projectDir}/assets/stub/domains.tsv domains.tsv
    """
}

process STRUCTURE_HITS {
    tag   "${sample}"
    label 'cpu'
    conda "${params.foldseek_env}"
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'structure_hits.tsv'
    input:
    tuple val(sample), path(struct_tsv), path(struct_dir)
    output:
    tuple val(sample), path("structure_hits.tsv"), emit: feat
    script:
    """
    # foldseek easy-search the folded structures against the CAZyme3D reference DB
    foldseek easy-search ${struct_dir} ${params.cazyme3d_db} aln.tsv tmp_fs \\
        --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,lddt' \\
        -e 1e-3 || : > aln.tsv
    feature_converters.py structure_hits --aln aln.tsv --out structure_hits.tsv --topk ${params.struct_topk}
    """
    stub:
    """
    cp ${projectDir}/assets/stub/structure_hits.tsv structure_hits.tsv
    """
}

process LOCALIZATION {
    tag   "${sample}"
    label 'cpu'
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'localization.tsv'
    input:
    tuple val(sample), path(signalp_tsv)
    output:
    tuple val(sample), path("localization.tsv"), emit: feat
    script:
    """
    feature_converters.py localization --signalp ${signalp_tsv} --out localization.tsv
    """
    stub:
    """
    cp ${projectDir}/assets/stub/localization.tsv localization.tsv
    """
}

process PHYSICOCHEM {
    tag   "${sample}"
    label 'cpu'
    conda "${params.dbcan_env}"
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'physicochem.tsv'
    input:
    tuple val(sample), path(faa)
    output:
    tuple val(sample), path("physicochem.tsv"), emit: feat
    script:
    """
    feature_converters.py physicochem --faa ${faa} --out physicochem.tsv
    """
    stub:
    """
    cp ${projectDir}/assets/stub/physicochem.tsv physicochem.tsv
    """
}

process CLEAN_EC {
    tag   "${sample}"
    label 'gpu'
    conda "${params.clean_env}"
    publishDir { "${params.outdir}/cazyme_advanced/features/${sample}" }, mode: 'copy', pattern: 'ec_prediction.tsv'
    input:
    tuple val(sample), path(faa)
    output:
    tuple val(sample), path("ec_prediction.tsv"), emit: feat
    script:
    """
    # CLEAN inference (maxsep); wrapper writes results/<name>_maxsep.csv
    run_clean.sh ${faa} clean_out.csv || : > clean_out.csv
    feature_converters.py ec_prediction --clean clean_out.csv --out ec_prediction.tsv
    """
    stub:
    """
    cp ${projectDir}/assets/stub/ec_prediction.tsv ec_prediction.tsv
    """
}
