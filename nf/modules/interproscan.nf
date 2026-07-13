// InterProScan tier — produces the funcscan-layout InterProScan TSV that
// BioForge's baseline ingester (bioforge-ingest) discovers, populating the
// Gene Ontology card and the InterPro-domains table on the gene deep-dive page.
//
// Two paths, chosen at run time (no config edit needed to switch):
//   * REAL:     if params.interproscan_sh points at an interproscan.sh install,
//               run it natively (-f tsv -goterms -pa) for full member-DB + GO coverage.
//   * FALLBACK: otherwise derive Analysis=Pfam signatures + GO terms from the Pfam
//               domains the pipeline already computed (domains.tsv), joined against
//               the vendored pfam2go map — offline, no tens-of-GB install. This keeps
//               the "git clone and run" promise while still populating GO/InterPro.
//
// Output lands at funcscan/protein_annotation/interproscan/<sample>_interproscan_faa.tsv,
// exactly where bioforge.ingest.discover looks for it.

process INTERPROSCAN {
    tag   "${sample}"
    label 'cpu'
    conda "${params.dbcan_env}"
    publishDir "${params.outdir}/funcscan", mode: 'copy',
        saveAs: { fn -> fn.startsWith('funcscan/') ? fn.substring('funcscan/'.length()) : fn }

    input:
    tuple val(sample), path(faa), path(domains_tsv)

    output:
    tuple val(sample), path("funcscan/protein_annotation/interproscan/${sample}_interproscan_faa.tsv"), emit: tsv

    script:
    def ips   = params.interproscan_sh ?: ''
    def p2ipr = params.pfam2interpro ? "--pfam2interpro ${params.pfam2interpro}" : ''
    """
    mkdir -p funcscan/protein_annotation/interproscan
    OUT=funcscan/protein_annotation/interproscan/${sample}_interproscan_faa.tsv
    if [ -n "${ips}" ] && [ -x "${ips}" ]; then
        echo "[interproscan] real InterProScan: ${ips}"
        "${ips}" -i ${faa} -f tsv -goterms -pa -o "\$OUT" -T ips_tmp || {
            echo "[interproscan] real run failed -> Pfam->GO fallback" >&2
            pfam_to_interproscan.py --domains-tsv ${domains_tsv} \\
                --pfam2go ${params.pfam2go} ${p2ipr} --faa ${faa} --out "\$OUT"
        }
    else
        echo "[interproscan] interproscan.sh not configured -> deriving GO from Pfam domains"
        pfam_to_interproscan.py --domains-tsv ${domains_tsv} \\
            --pfam2go ${params.pfam2go} ${p2ipr} --faa ${faa} --out "\$OUT"
    fi
    """

    stub:
    """
    mkdir -p funcscan/protein_annotation/interproscan
    cp ${projectDir}/assets/stub/interproscan.tsv \\
       funcscan/protein_annotation/interproscan/${sample}_interproscan_faa.tsv
    """
}
