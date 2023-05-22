#!/usr/bin/env nextflow
/*
Zymo smrnaseq pipeline, adapted from nf-core/smrnaseq Analysis Pipeline (https://github.com/nf-core/smrnaseq).
Many changes were made including migration to DSL 2
*/

// DSL 2
nextflow.enable.dsl=2

// Load functions
include { help_message } from ('./libs/help_message')
include { parse_presets } from ('./libs/parse_presets')
include { setup_channel } from ('./libs/setup_channel')
include { collect_summary } from ('./libs/collect_summary')

// Show help message
if (params.help){
    help_message()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
// Genome options
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}
params.hairpin = params.genome ? params.genomes[ params.genome ].hairpin ?: false : false
params.ncRNA_fa = params.genome ? params.genomes[params.genome].ncRNA_fa ?: false : false
params.ncRNA_fa_bowtie = params.genome ? params.genomes[ params.genome ].ncRNA_fa_bowtie ?: false : false
params.idmaptoRNAtype = params.genome ? params.genomes[ params.genome ].idmaptoRNAtype ?: false : false
params.tx2gene = params.genome ? params.genomes[ params.genome ].tx2gene ?: false : false
params.mirtrace_species = params.genome ? params.genomes[ params.genome ].mirtrace_species ?: false : false
if( !params.mirtrace_species && params.genome != "artificial" ) {
    exit 1, "Reference species for miRTrace is not defined."
}
if (params.mirna_gtf) {
    params.mirna_gtf_path = params.mirna_gtf
} else if (params.genome && params.genomes[ params.genome ].mirna_gtf) {
    params.mirna_gtf_path = params.genomes[ params.genome ].mirna_gtf
} else if (params.mirtrace_species) {
    params.mirna_gtf_path = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${params.mirtrace_species}.gff3"
} else {
    params.mirna_gtf_path = false
}

// AWSBatch sanity checking
if( workflow.profile == 'awsbatch') {
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Parse protocol presets
(clip_R1, three_prime_clip_R1, adapter, protocol) = parse_presets(params.protocol)

// --outdir Remove potential trailing slash at the end
outdir = params.outdir - ~/\/*$/

// Maximum and minimum length check
if (params.max_length && params.max_length < params.min_length) {
    exit 1, "Maximum length has been set smaller than minimum read length."
}

/*
 * SET & VALIDATE INPUT CHANNELS
 */
hairpin = setup_channel(params.hairpin, "hairpin FASTA file", false, "alignment against hairpin and following steps will be skipped.")
ncRNA_fa = setup_channel(params.ncRNA_fa, "RNA Types FASTA file", false, "Quantification of RNA Types and following steps will be skipped.")
ncRNA_fa_bowtie = setup_channel(params.ncRNA_fa_bowtie ? "${params.ncRNA_fa_bowtie}*" : false,  "RNA Types Bowtie Index", false, "Alignment and quantification of RNA Types and following steps will be skipped.")
idmaptoRNAtype = setup_channel(params.idmaptoRNAtype ? "${params.idmaptoRNAtype}" : false, "RNA Types to IDs Map", false, "Quantification of RNA Types and following steps will be skipped.")
tx2gene = setup_channel(params.tx2gene ? "${params.tx2gene}" : false, "RNA Transcript IDs to Genes Map", false, "Differential Expression of RNA Types and following steps will be skipped.")
mirna_gtf = setup_channel(params.mirna_gtf_path, "miRNA GTF/GFF file", false, "mirtop step will be skipped.")
design = setup_channel(params.design, "design CSV file", true, "")
comparisons = setup_channel(params.comparisons, "comparison CSV file", false, "all pairwise comparisons will be carried out.")
summary_header = Channel.fromPath("$baseDir/assets/workflow_summary_header.txt")
fastqc_limits = Channel.fromPath("$baseDir/assets/FastQC_limits.txt")
multiqc_config = Channel.fromPath("$baseDir/assets/multiqc_config.yaml")

/*
 * COLLECT SUMMARY & LOG
 */
log.info "Aladdin smrnaseq v${workflow.manifest.version}"
def summary = collect_summary(params, workflow)
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
// Create a channel for settings to show in MultiQC report
info_to_report = ['Genome', 'Protocol', 'Min Trimmed Length', 'Max Trimmed Length', "3' adapter", "Trim 5' R1", "Trim 3' R1", 'isomiRs FDR cutoff', 'isomiRs Log2FC cutoff']
summary_to_report = summary.findAll { it.key in info_to_report }
Channel.from( summary_to_report.collect{ [it.key, it.value] } )
       .map { k,v -> "            <dt>$k</dt><dd><samp>$v</samp></dd>" }
       .collectFile(name: 'summary.txt', newLine: true, sort: 'index')
       .set { workflow_summary_to_report }
// Save workflow summary plain text
Channel.from( summary.collect{ [it.key, it.value] } )
       .map { k,v -> "${k.padRight(21)} : $v" }
       .collectFile(name: "${outdir}/pipeline_info/workflow_summary.txt", newLine: true, sort: 'index')

/*
 * PROCESS DEFINITION
 */
include { check_design } from "./processes/check_design"
include { make_bowtie_index } from "./processes/make_bowtie_index" addParams(
    publish_dir: "${outdir}/bowtie/reference",
    mirtrace_species: params.mirtrace_species,
    save_reference: params.save_reference
)
include { fastqc } from "./processes/fastqc" addParams(
    publish_dir: "${outdir}/fastqc",
    skip_fastqc: params.skip_fastqc
)
include { trim_galore } from "./processes/trim_galore" addParams(
    outdir: outdir,
    clip_R1: clip_R1,
    three_prime_clip_R1: three_prime_clip_R1,
    adapter: adapter,
    protocol: protocol,
    min_length: params.min_length,
    max_length: params.max_length,
    trim_nextseq: params.trim_nextseq,
    save_trimmed: params.save_trimmed
)
include { mirtrace } from "./processes/mirtrace" addParams(
    publish_dir: "${outdir}/miRTrace",
    mirtrace_species: params.mirtrace_species,
    protocol: protocol
)
include { collapse_reads } from "./processes/collapse_reads" addParams(
    genome: params.genome
)
include { bowtie as bowtie_hairpin } from "./processes/bowtie" addParams(
    publish_dir: "${outdir}/bowtie/miRBase_collapsed",
    align_mode: "-a"
)
include { bowtie as bowtie_rnacounts } from "./processes/bowtie" addParams(
    publish_dir: "${outdir}/bowtie/RNAtypecounts",
    align_mode: "-a"
)
include { make_rsem_ref } from "./processes/make_rsem_ref" addParams(
   publish_dir: "${outdir}/rsem/reference",
)
include { rsem } from "./processes/rsem" addParams(
   publish_dir: "${outdir}/rsem",
)
include { tximport_deseq } from "./processes/tximport_deseq" addParams(
   publish_dir: "${outdir}/tximport_deseq"
)
include { quantify_rnacounts } from "./processes/quantify_rnacounts" addParams(
   publish_dir: "${outdir}/quantify_rnacounts",
   max_length: params.max_length
)
include { mirtop } from "./processes/mirtop" addParams(
    outdir: outdir,
    mirtrace_species: params.mirtrace_species
)
include { isomirs } from "./processes/isomirs" addParams(
    publish_dir: "${outdir}/isomiRs",
    isomirs_fdr: params.isomirs_fdr,
    isomirs_lfc: params.isomirs_lfc
)
include { software_versions } from "./processes/software_versions" addParams(
    publish_dir: "${outdir}/pipeline_info",
    pipeline_version: workflow.manifest.version,
    nextflow_version: workflow.nextflow.version
)
include { multiqc } from "./processes/multiqc" addParams(
    publish_dir: "${outdir}/MultiQC",
    skip_multiqc: params.skip_multiqc,
    run_name : summary["Run Name"],
    isomirs_fdr: params.isomirs_fdr
)
include { summarize_downloads } from "./processes/summarize_downloads" addParams(
    publish_dir: "${outdir}/download_data"
)

/*
 * WORKFLOW DEFINITION
 */
workflow {
    check_design(design, comparisons.ifEmpty([]))
    check_design.out.checked_design
        .splitCsv( header: true )
        .map { row -> [ row["sample"], [ file(row["read_1"]) ] ] }
        .set { reads }
    make_bowtie_index(hairpin)
    fastqc(reads, fastqc_limits.collect())
    trim_galore(reads)
    mirtrace(trim_galore.out.reads.collect())
    collapse_reads(trim_galore.out.reads)
    bowtie_hairpin(collapse_reads.out.fastq, make_bowtie_index.out.index.collect())
    bowtie_rnacounts(trim_galore.out.reads, ncRNA_fa_bowtie.collect())
    make_rsem_ref(ncRNA_fa.collect())
    rsem(bowtie_rnacounts.out.bam, make_rsem_ref.out.rsem_ref.collect())
    tximport_deseq(rsem.out.rsem_quant.collect(), check_design.out.checked_design, tx2gene, comparisons.ifEmpty([]))
    quantify_rnacounts(bowtie_rnacounts.out.bam.collect(), tximport_deseq.out.tximporttsv, idmaptoRNAtype, rsem.out.rsem_logs.collect())
    mirtop(bowtie_hairpin.out.bam.collect(), make_bowtie_index.out.fasta, mirna_gtf)
    isomirs(mirtop.out.isomir_tsv, check_design.out.checked_design, comparisons.ifEmpty([]))
    software_versions()
    multiqc(multiqc_config, \
            fastqc.out.report.collect().ifEmpty([]), \
            trim_galore.out.report.collect(), \
            mirtrace.out.report.mix(mirtrace.out.general_stats).collect().ifEmpty([]), \
            quantify_rnacounts.out.miRNA_perc_gs.mix(quantify_rnacounts.out.RNAtypechart).collect().ifEmpty([]), \
            isomirs.out.report.collect().ifEmpty([]), \
            tximport_deseq.out.report.collect().ifEmpty([]), \
            software_versions.out.report.collect(), \
            summary_header, \
            workflow_summary_to_report)
    isomirs_locations =  isomirs.out.download.flatten().map{ "${outdir}/isomiRs/" + it.getName() }
    deseq_locations = tximport_deseq.out.download.flatten().map{ "${outdir}/tximport_deseq/" + it.getName()}
    report_locations = multiqc.out.report.map{ "${outdir}/MultiQC/" + it.getName() }
    aladdin_locations = isomirs_locations.mix(deseq_locations, report_locations).collectFile(name: "${outdir}/download_data/file_locations.txt", newLine: true)
    summarize_downloads(aladdin_locations, check_design.out.checked_design)
}

/*
 * LOG ON COMPLETION
 */
workflow.onComplete {
    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "Warning, pipeline completed, but with errored process(es)"
      log.info "Number of ignored errored process(es) : ${workflow.stats.ignoredCount}"
      log.info "Number of successfully ran process(es) : ${workflow.stats.succeedCount}"
    }
    if(workflow.success){
        log.info "[smrnaseq] Pipeline completed successfully"
    } else {
        log.info "[smrnaseq] Pipeline completed with errors"
    }
}
