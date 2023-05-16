// Sample comparison and statistical analysis with isomiRs
params.publish_dir = "isomiRs"
params.isomirs_fdr = 0.05
params.isomirs_lfc = 0.585

process isomirs {
    label 'process_low'
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
    path mirtop
    path samplegrouping
    path comparisons

    output:
    path "*.{xlsx,jpg}", emit: download
    path "*.tsv", emit: report optional true

    script:
    """
    isomiRs.r $mirtop $samplegrouping $params.isomirs_fdr $params.isomirs_lfc $comparisons 
    """
}
