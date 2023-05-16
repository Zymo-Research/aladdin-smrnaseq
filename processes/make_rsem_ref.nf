params.publish_dir="rsem/reference"

process make_rsem_ref{
     publishDir "$params.publish_dir", mode:'copy'

     input:
     path ref_fa

     output:
     path 'hsa_rsemref*', emit: rsem_ref

     script:
     """
     rsem-prepare-reference -p $task.cpus ${ref_fa} "hsa_rsemref"
     """
}
