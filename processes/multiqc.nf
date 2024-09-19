// Generate report with MultiQC
params.publish_dir = "MultiQC"
params.skip_multiqc = false
params.run_name = false
params.isomirs_fdr = 0.05
params.ensembl_link = false
params.mirbase_link = false

process multiqc {
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    path ('original_multiqc_config.yaml')
    path ('fastqc/*')
    path ('trim_galore/*')
    path ('mirtrace/*')
    path ('quantify_rnacounts/*')
    path ('isomirs/*')
    path ('tximport/*')
    path ('software_versions/*')
    path summary_header
    path workflow_summary

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"
    //file "multiqc_plots"

    script:
    rtitle = params.run_name ? "--title \"${params.run_name}\"" : ''
    rfilename = params.run_name ? "--filename " + params.run_name.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    ensembl_link = params.ensembl_web ? "ensembl_link_prefix: ${params.ensembl_web}" : ''
    mirbase_link = params.mirbase_web ? "mirbase_link_prefix: ${params.mirbase_web}" : ''
    """
    cat $summary_header $workflow_summary > workflow_summary_mqc.yaml
    echo '    </dl>' >> workflow_summary_mqc.yaml
    rm $summary_header $workflow_summary
    cat original_multiqc_config.yaml > multiqc_config.yaml
    rm original_multiqc_config.yaml
    echo '$ensembl_link' >> multiqc_config.yaml
    echo '$mirbase_link' >> multiqc_config.yaml
    echo 'isomiRs_alpha: ${params.isomirs_fdr}' >> multiqc_config.yaml
    multiqc . -f $rtitle $rfilename --config multiqc_config.yaml -m mirtrace -m Trim_Galore \\
              -m fastqc -m custom_content -m isomiRs -m deseq2_rnatypes -m plot_sample_distance -m plot_gene_heatmap
    """
}
