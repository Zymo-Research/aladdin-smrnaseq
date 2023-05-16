params.publish_dir = "downsample"
params.downsample_num = false

process downsample {
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
    tuple val(name), path(reads)

    when:
    params.downsample_num 

    output:
    tuple val(name), path('*_R1.fastq.gz'), emit: ds_reads

    script:
    """
    readnum=\$((\$(zcat $reads | wc -l) / 4))
    if ((\$readnum > $params.downsample_num))
    then
    seqtk sample -s1000 $reads $params.downsample_num > ${name}_downsample_R1.fastq
    gzip ${name}_downsample_R1.fastq
    else
    [ ! -f ${name}_R1.fastq.gz ] && ln -s $reads ${name}_R1.fastq.gz
    fi
    """
}
