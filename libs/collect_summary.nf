def collect_summary(params, workflow) {
    def summary = [:]
    if (workflow.revision) {
        summary['Pipeline Release'] = workflow.revision
    }          
    // This has the bonus effect of catching both -name and --name
    def run_name = params.name
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        run_name = workflow.runName
    }
    
    presets = [ "Zymo-Seq_miRNA" : [ "clip_R1":1, "three_prime_clip_R1":0, "adapter":"TGGAATTCTCGGGTGCCAAGG", "protocol":"illumina" ],
                "illumina"       : [ "clip_R1":0, "three_prime_clip_R1":0, "adapter":"TGGAATTCTCGGGTGCCAAGG", "protocol":"illumina" ],
                "nextflex"       : [ "clip_R1":4, "three_prime_clip_R1":4, "adapter":"TGGAATTCTCGGGTGCCAAGG", "protocol":"nextflex" ],
                "qiaseq"         : [ "clip_R1":0, "three_prime_clip_R1":0, "adapter":"AACTGTAGGCACCATCAAT", "protocol":"qiaseq" ],
                "cats"           : [ "clip_R1":3, "three_prime_clip_R1":0, "adapter":"AAAAAAAA", "protocol":"cats" ]
    ]

    summary['Run Name']              = run_name ?: workflow.runName
    summary['Design']                = params.design
    summary['Genome']                = params.genome
    summary['Protocol']              = params.protocol
    summary['Min Trimmed Length']    = params.min_length
    summary['Max Trimmed Length']    = !params.max_length ? "None" : params.max_length
    summary['miRBase hairpin']       = params.hairpin
    summary['Save Reference']        = params.save_reference ? 'Yes' : 'No'
    summary['miRTrace species']      = params.mirtrace_species
    summary['Group Comparisons']     = params.comparisons
    summary['isomiRs FDR cutoff']    = params.isomirs_fdr
    summary['isomiRs Log2FC cutoff'] = params.isomirs_lfc
    summary['Output dir']            = params.outdir
    summary['Launch dir']            = workflow.launchDir
    summary['Working dir']           = workflow.profile == 'awsbatch' ? "s3:/${workflow.workDir}" : workflow.workDir
    summary['Config Profile']        = workflow.profile
    summary['Max Resources']         = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
    if (workflow.profile == 'awsbatch') {
        summary['AWS Region']        = params.awsregion
        summary['AWS Queue']         = params.awsqueue
    }
    summary["Trim 5' R1"]            = presets[params.protocol]["clip_R1"]
    summary["Trim 3' R1"]            = presets[params.protocol]["three_prime_clip_R1"]
    summary["3' adapter"]            = presets[params.protocol]["adapter"]
    return summary
}
