// Baseline dbCAN tier — run_dbcan V5 (HMMER + dbCAN-sub + DIAMOND) on the input
// protein FASTA, reshaped into the nf-core/funcscan layout that BioForge's
// baseline ingester (bioforge-ingest) discovers. This makes the pipeline
// self-contained: one protein FASTA in -> both baseline and advanced calls out,
// with no separate funcscan run required.
//
// Real execution uses the project venv on met (run_dbcan installed there) + the
// current dbCAN DB. The stub reproduces the published funcscan layout with a tiny
// canned overview so the DAG + downstream ingest can be proven with no tools.

process BASELINE_DBCAN {
    tag   "${sample}"
    label 'cpu'
    conda "${params.dbcan_env}"
    publishDir "${params.outdir}/funcscan", mode: 'copy',
        saveAs: { fn -> fn.startsWith('funcscan/') ? fn.substring('funcscan/'.length()) : fn }

    input:
    tuple val(sample), path(faa)

    output:
    tuple val(sample), path("funcscan/**"), emit: tree
    tuple val(sample), path("funcscan/cazyme/dbcan/cazyme_annotation/${sample}_overview.tsv"), emit: overview

    script:
    """
    mkdir -p rundbcan
    run_dbcan CAZyme_annotation \\
        --mode protein \\
        --input_raw_data ${faa} \\
        --db_dir ${params.dbcan_db} \\
        --output_dir rundbcan \\
        --methods diamond,hmm,dbCANsub
    emit_baseline_funcscan.py \\
        --overview rundbcan/overview.tsv \\
        --faa ${faa} \\
        --sample ${sample} \\
        --outdir funcscan
    """

    stub:
    """
    mkdir -p funcscan
    # Canned overview covering ALL protein ids the advanced/feature stubs use
    # (demo_p01..demo_p10) so a stub run ingests a fully self-consistent
    # baseline+advanced demo with no "unmatched to a gene" warnings.
    printf 'Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n' > overview_stub.tsv
    printf 'demo_p01\t3.2.1.-\tGH43_26(20-350)\tGH43_26\tGH43_26\t3\tGH43_26\t-\n' >> overview_stub.tsv
    printf 'demo_p02\t-\t-\t-\tGT4\t1\tGT4\t-\n' >> overview_stub.tsv
    printf 'demo_p03\t3.2.1.-\tGH31_1(30-480)\tGH31_1\tGH31_1\t3\tGH31_1\t-\n' >> overview_stub.tsv
    printf 'demo_p04\t3.2.1.4\tGH5(40-410)\tGH5_5\tGH5\t3\tGH5\t-\n' >> overview_stub.tsv
    printf 'demo_p05\t-\tGH13(15-460)\tGH13_1\tGH13\t3\tGH13\t-\n' >> overview_stub.tsv
    printf 'demo_p06\t2.4.1.-\t-\t-\tGT2\t1\tGT2\t-\n' >> overview_stub.tsv
    printf 'demo_p07\t3.2.1.8\tGH10(25-330)\tGH10\tGH10\t3\tGH10\t-\n' >> overview_stub.tsv
    printf 'demo_p08\t3.2.1.21\tGH3(30-720)\tGH3\tGH3\t3\tGH3\t-\n' >> overview_stub.tsv
    printf 'demo_p09\t-\tCBM1(1-36)\t-\tCBM1\t2\tCBM1\t-\n' >> overview_stub.tsv
    printf 'demo_p10\t3.2.1.91\tGH6(28-390)\tGH6\tGH6\t3\tGH6\t-\n' >> overview_stub.tsv
    {
      for p in 01 02 03 04 05 06 07 08 09 10; do
        printf '>demo_p%s\\nMKTAYIAKQRVVPASGTNETDDAPMLLRIVIILAALARVSAGGMYLLNSVMNISWLQDVDIDLR\\n' "\$p"
      done
    } > stub.faa
    emit_baseline_funcscan.py --overview overview_stub.tsv --faa stub.faa \\
        --sample ${sample} --outdir funcscan
    """
}
