// FastQC before trimming
params.publish_dir = "fastqc"
params.skip_fastqc = false

process fastqc {
    label "process_low"
    tag "$name"
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    !params.skip_fastqc

    input:
    tuple val(name), path(reads)
    path limits

    output:
    path '*_fastqc.{zip,html}', emit: report

    script:
    """
    [ ! -f  ${name}_R1.fastq.gz ] && ln -s $reads ${name}_R1.fastq.gz
    fastqc --quiet --threads $task.cpus -l $limits ${name}_R1.fastq.gz
    """
}