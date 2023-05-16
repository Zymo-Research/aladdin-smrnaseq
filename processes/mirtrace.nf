// miRNA QC with miRTrace
params.mirtrace_species = "hsa"
params.publish_dir = "miRTrace"
params.protocol = "illumina"

process mirtrace {
    label 'process_low'    
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    params.mirtrace_species
      
    input:
    path reads

    output:
    path 'mirtrace/mirtrace*.{tsv,json,html}', emit: report
    path 'mirtrace/qc_passed_reads*'
    path '*mirtrace_mqc.json', emit: general_stats

    script:
    """
    for i in $reads
    do
        path=\$(realpath \${i})
        prefix=\$(echo \${i} | sed -e "s/_trimmed.fq.gz//")
        echo \$path","\$prefix
    done > mirtrace_config

    mirtrace qc \\
        --species $params.mirtrace_species \\
        --protocol $params.protocol \\
        --config mirtrace_config \\
        --write-fasta \\
        --output-dir mirtrace \\
        -t $task.cpus \\
        --force

    miRNA_percent_mqc.py mirtrace/mirtrace-stats-rnatype.tsv
    miRNA_species_mqc.py mirtrace/mirtrace-stats-mirna-complexity.tsv
    """
 }