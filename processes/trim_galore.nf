// Trimming
params.outdir = "./"
params.save_trimmed = false
params.clip_R1 = 0
params.three_prime_clip_R1 = 0
params.adapter = "TGGAATTCTCGGGTGCCAAGG"
params.protocol = "illumina"
params.min_length = 18
params.max_length = false
params.trim_nextseq = 0


process trim_galore {
    label 'process_medium'
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (params.save_trimmed && filename.endsWith("fq.gz")) filename
            else null
        } 

    input:
    tuple val(name), path(reads)

    output:
    path '*fq.gz', emit: reads
    path '*trimming_report.txt', emit: report
    path "*_fastqc.{zip,html}", emit: fastqc
    path "${name}_R1.fastq.gz", includeInputs: true, emit: download

    script:
    tg_min_length = "--length ${params.min_length}"
    tg_max_length = !params.max_length ? "": "--max_length ${params.max_length}"
    c_r1 = params.clip_R1 > 0 ? "--clip_R1 ${params.clip_R1}" : ''
    tpc_r1 = params.three_prime_clip_R1 > 0 ? "--three_prime_clip_R1 ${params.three_prime_clip_R1}" : ''
    tpa = (params.protocol == "qiaseq" | params.protocol == "cats") ? "--adapter ${params.adapter}" : '--small_rna'
    nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
    """
    [ ! -f  ${name}_R1.fastq.gz ] && ln -s $reads ${name}_R1.fastq.gz
    trim_galore $tpa $tg_min_length $tg_max_length $c_r1 $tpc_r1 $nextseq -j $task.cpus --gzip --fastqc --basename $name ${name}_R1.fastq.gz
    """
}
