// Align reads to reference with bowtie
params.align_mode = "-a"
params.publish_dir = "bowtie"

process bowtie {
    label 'process_medium'
    tag "$prefix"
    publishDir "${params.publish_dir}", mode: 'copy'
    
    input:
    path reads
    path index

    output:
    path '*.bam', emit: bam

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() -  ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    read_cmd = reads.toString().endsWith('.gz') ? 'zcat' : 'cat'
    """
    bowtie \\
        $index_base \\
        -p ${task.cpus} \\
        -t \\
        -m 50 \\
        ${params.align_mode} \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        -q <($read_cmd $reads) \\
        -S \\
        | samtools view -bS - > ${prefix}.bam
    """
}