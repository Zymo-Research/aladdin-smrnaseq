params.publish_dir = "quantify_rnacounts"
params.max_length = false

process quantify_rnacounts {
    
    publishDir "${params.publish_dir}", mode: 'copy'
 
    input:
    path bamfiles
    path tximportdata
    path idmaptornatype
    path rsemlogs

    output:
    path "*.csv", emit: RNAcountchart
    path "perc_miRNA_quantifyrna_mqc.json", emit: miRNA_perc_gs
    path "RNAtypechart_mqc.json", emit: RNAtypechart

    script:
    bamlist = bamfiles.toString()
    qr_max_length = !params.max_length ? "" : "${params.max_length}"
    """
    quantify_rnacounts.py -t "${tximportdata}" -i "${idmaptornatype}" -f "${bamlist} -r ${rsemlogs}"
    quantifyrna_miRNAperc_mqc.py -d "RNAcounts.csv" 
    displayRNATypeCounts_mqc.py -d "RNAcounts.csv" -m "$qr_max_length"
    """
}








