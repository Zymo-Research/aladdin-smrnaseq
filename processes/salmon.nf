params.publish_dir="salmon"
params.deliver_auxfiles = false

process salmon {
        label 'process_medium'
        if (params.deliver_auxfiles) {

            publishDir "$params.publish_dir/aux_info", mode: 'copy',
               saveAs: {filename ->
                   if (filename.indexOf("aux_info") >0) "$filename"
                   else null
            }
        }

	publishDir "$params.publish_dir", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("quant.log") > 0 ) "logs/$filename"
                else if (filename.indexOf("quant.sf") > 0 ) "quant/$filename"
                else if (filename.indexOf("_info.json") > 0 ) "logs/run_info/$filename"
                else null
            }

	input:
	path bams
	path ncRNA_fa

	output:
	path '*quant.sf', emit: salmonquant
        path '*quant.log', emit: salmonlog
        path '*_info.json'
        path '*aux_info/*'

	script:
	prefix = bams.toString() - ".bam"	
	"""
	salmon quant -t $ncRNA_fa -l A -a $bams -p $task.cpus -o salmonquant
        mv salmonquant/quant.sf ${prefix}_quant.sf
        mv salmonquant/logs/salmon_quant.log ${prefix}_salmon_quant.log
        mv salmonquant/cmd_info.json ${prefix}_cmd_info.json
        mv salmonquant/aux_info/meta_info.json ${prefix}_meta_info.json
        mv salmonquant/aux_info/ ${prefix}_aux_info/
	"""
}

