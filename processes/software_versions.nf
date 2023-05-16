// Gather software versions
params.publish_dir = "pipeline_info"
params.pipeline_version = "1.0.0"
params.nextflow_version = "20.07.1"

process software_versions {
   label 'process_low'
   publishDir "${params.publish_dir}", mode: 'copy',
   saveAs: { it == "software_versions.csv" ? it : null }
   
   output:
   path 'software_versions_mqc.yaml', emit: report
   path 'software_versions.csv'

   script:
   """
   echo $params.pipeline_version > v_pipeline.txt
   echo $params.nextflow_version > v_nextflow.txt
   echo \$(R --version 2>&1) > v_R.txt
   fastqc --version > v_fastqc.txt
   trim_galore --version > v_trim_galore.txt
   bowtie --version > v_bowtie.txt
   samtools --version > v_samtools.txt
   seqcluster --version > v_seqcluster.txt
   fasta_formatter -h > v_fastx.txt
   mirtop --version > v_mirtop.txt
   mirtrace --version > v_mirtrace.txt
   rsem-calculate-expression --version > v_rsem.txt
   Rscript -e 'write(x=as.character(packageVersion("tximport")), file="v_tximport.txt")'
   Rscript -e 'write(x=as.character(packageVersion("isomiRs")), file="v_isomiRs.txt")'
   Rscript -e 'write(x=as.character(packageVersion("DESeq2")), file="v_DESeq2.txt")'
   scrape_software_versions.py > software_versions_mqc.yaml
   """
}
