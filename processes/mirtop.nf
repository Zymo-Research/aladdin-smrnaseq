// Annotate and tally miRNA and isomiRs with mirtop
params.mirtrace_species = "hsa"
params.outdir = "./"

process mirtop {

    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path input
    path hairpin
    path gtf
    
    output:
    path "mirtop/mirtop.gff"
    path "mirtop/mirtop.tsv"
    path "mirtop/mirna.tsv"
    path "mirtop/mirtop_rawData.tsv", emit: isomir_tsv
    path "mirtop/mirtop_stats.*"
    
    script:
    """
    mirtop gff --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species $input
    mirtop stats -o mirtop mirtop/mirtop.gff
    mirtop counts --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species --add-extra --gff mirtop/mirtop.gff
    mirtop export --format isomir --hairpin $hairpin --gtf $gtf --sps $params.mirtrace_species -o mirtop mirtop/mirtop.gff
    collapse_mirtop.r mirtop/mirtop.tsv
    """   
}
