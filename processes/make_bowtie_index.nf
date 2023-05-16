// Make bowtie index for hairpin fasta
params.mirtrace_species = "hsa"
params.save_reference = false
params.publish_dir = "bowtie/reference"

process make_bowtie_index {
    publishDir "${params.publish_dir}", mode: "copy",
               saveAs: { params.save_reference ? it : null }

    input:
    path hairpin

    output:
    path "hairpin_idx.*", emit: index
    path "hairpin_idx.fa", emit: fasta

    script:
    """
    seqkit grep -r --pattern \".*${params.mirtrace_species}-.*\" $hairpin > hairpin_sps.fa
    seqkit seq --rna2dna hairpin_sps.fa > hairpin_igenome.fa
    fasta_formatter -w 0 -i hairpin_igenome.fa -o hairpin_idx.fa
    bowtie-build hairpin_idx.fa hairpin_idx
    """
}