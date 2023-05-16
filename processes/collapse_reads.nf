// Collapse reads with same sequences
params.genome = false

process collapse_reads {
    tag "$prefix"

    when:
    !(params.genome == "artificial")
    
    input:
    path reads

    output:
    path 'collapsed/*.fastq', emit: fastq

    script:
    prefix = reads.toString() - '_trimmed.fq.gz'
    """
    seqcluster collapse -f $reads -m 1 --min_size 15 -o collapsed
    mv collapsed/${prefix}_trimmed_trimmed.fastq collapsed/${prefix}.fastq
    """
}
