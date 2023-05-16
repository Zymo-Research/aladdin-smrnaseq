params.publish_dir="rsem"

process rsem{
     label 'process_medium'
     publishDir "$params.publish_dir", mode:'copy'

     input:
     path bams
     path rsem_ref

     output:
     path '*.isoforms.results', emit: rsem_quant
     path '*rsem.log' , emit: rsem_logs

     script:
     prefix = bams.toString() - ".bam"
     """
     rsem-calculate-expression -p $task.cpus --no-bam-output --seed-length 18  --alignments --bam $bams "hsa_rsemref" ${prefix}
     mv .command.out ${prefix}_rsem.log
     """
}
