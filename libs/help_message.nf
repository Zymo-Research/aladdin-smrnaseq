def help_message() {
    log.info "Zymo smrnaseq v${workflow.manifest.version}"
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --design "sample_sheet.csv" --genome GRCh38 -profile docker

    Mandatory arguments:
      --design                      Path to a CSV file with sample labels, sample groupings and input FASTQ file locations
      --genome                      Name of reference genome. See conf/igenomes.conf for list of supported genomes. Artificial polymers are also supported with option "artificial", which skips alignments and subsequent steps.
      --protocol                    Library preparation protocol. Can be set as "illumina", "nextflex", "qiaseq", "cats", or "Zymo-Seq_miRNA"

    References
      --save_reference              Save the generated reference files the the Results directory. Default: false
      --mirna_gtf                   GFF/GTF file with coordinates positions of precursor and miRNAs. Pre-defined when '--genome' is specified. Use this to override that.

    Trimming options
      --min_length [int]            Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18
      --max_length [int]            Discard reads that are longer than [int] after adapter trimming. Default: false
      --save_trimmed                Save trimmed FASTQ file intermediates

    QC
      --skip_fastqc                 Skip FastQC
      --skip_multiqc                Skip MultiQC

    Comparisons
      --comparisons                 Path to a CSV file stating the sample groups you want to compare. If not provided, all pairwise comparisons will be carried out.
      --isomirs_fdr [float]         FDR cutoff for isomiRs, default is 0.05
      --isomirs_lfc [float]         Log2FoldChange for isomiRs, please use log2 of desired fold change cutoff, default is 0.585, which is log2(1.5)

    Other options
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on

    """.stripIndent()
}
