params.publish_dir = "tximport_deseq"
params.isomirs_fdr = 0.05
params.isomirs_lfc = 0.585

process tximport_deseq {
   label 'process_low'
   publishDir "${params.publish_dir}", mode : 'copy'

   input: 
   path rsem
   path samplegrouping
   path tx2gene
   path comparisons

   output:
   path 'tximport_filemap.tsv'
   path 'rsem_estimated_transcriptcounts.tsv', emit: tximporttsv
   path '*.tsv', emit: report
   path '*.{xlsx,jpg}', emit: download

   script:
   """
   for i in $rsem
   do
       path=\$(realpath \${i})
       echo \$path
   done > tximport_filemap.tsv

   tximport_deseq.r tximport_filemap.tsv $samplegrouping $params.isomirs_fdr $params.isomirs_lfc $tx2gene $comparisons
   """
} 
