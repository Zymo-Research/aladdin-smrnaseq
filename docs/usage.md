# nf-core/smrnaseq: Usage

## Table of contents
<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main Arguments](#main-arguments)
  * [`-profile`](#-profile)
  * [`--design`](#--design)
  * [`--protocol`](#--protocol)
* [Reference genomes](#reference-genomes)
  * [`--genome`](#--genome)
  * [`--mirna_gtf`](#--mirna_gtf)
  * [`--save_reference`](#--save_reference)
* [Trimming](#trimming)
  * [`--min_length`](#--min_length)
  * [`--max_length`](#--max_length)
  * [`--save_trimmed`](#--save_trimmed)
  * [`--trim_nextseq`](#--trim_nextseq)
* [Skipping QC steps](#skipping-qc-steps)
  * [`--skip_fastqc`](#--skip_fastqc)
  * [`--skip_multiqc`](#--skip_multiqc)
* [Comparisons](#comparisons)
  * [`--isomirs_fdr`](#--isomirs_fdr)
  * [`--isomirs_lfc`](#--isomirs_lfc)
  * [`--comparisons`](#--comparisons)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
  * [`-work-dir`](#-work-dir)
  * [`--outdir`](#--outdir)
  * [`--name`](#--name)
  * [`-resume`](#-resume)
<!-- TOC END -->


## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:

``` bash
nextflow run Zymo-Research/aladdin-smrnaseq \
    --design "<design CSV file>" \
    --genome GRCh38 \
    -profile awsbatch \
    -work-dir "<work dir on S3>" \
    --awsqueue "<SQS ARN>" \
    --outdir "<output dir on S3>"
```

This will launch the pipeline with the `awsbatch` configuration profile. See below for more information about profiles.

Because the repo is private, github credentials would be needed. You can set up that using [Nextflow SCM configuration file method](https://www.nextflow.io/docs/latest/sharing.html#scm-configuration-file), or you can simply `git clone` the repo to your computer/instance, and replace `Zymo-Research/aladdin-smrnaseq` with `main.nf` in the run command.

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Zymo-Research/aladdin-smrnaseq #if you have SCM config file set up
git pull Zymo-Research/aladdin-smrnaseq #if you prefer git
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Zymo-Research/aladdin-smrnaseq releases page](https://github.com/Zymo-Research/aladdin-smrnaseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main Arguments
### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

> While several profiles, such as `conda`, `docker`, and `singularity` can be used, only `awsbatch` has been extensively tested.
> We recommend you use `awsbatch` for most purposes.

### `--design`
Location of the design CSV file. The design CSV file must have the following format:

```
group,sample,read_1,read_2
Control,Sample1,s3://mybucket/this_is_s1_R1.fastq.gz,
Control,Sample2,s3://mybucket/this_is_s2_R1.fastq.gz,
Experiment,Sample3,s3://mybucket/that_is_s3_R1.fastq.gz,
Experiment,Sample4,s3://mybucket/that_be_s4_R1.fastq.gz,
```

Please note the following requirements for the design CSV file:
1. The header line must be present and cannot be changed.
2. Sample labels and group names must contain only alphanumerical characters or underscores, must start with a letter, and must not start with "R1" or "R2". Sample and group labels must also exclude phrases that will be automatically removed by MultiQC. Sample and group label terms unavailable for use in this pipeline can be found [here](https://github.com/ewels/MultiQC/blob/b936a7a6d7050f3edc1ceefe8ae6ecd93865bf66/multiqc/utils/config_defaults.yaml#L150-L284), in the MultiQC source code.
3. Only single-end sequencing is supported. Full S3 paths of R1 FASTQ files should be provided. R2 FASTQ files are ignored if provided.
4. If you do not wish to carry out statistical comparisons of samples, simply leave the group column blank (but keep the comma).

### `--protocol`
Protocol for constructing smRNA-seq libraries. Note that trimming parameters and 3' adapter sequence are pre-defined with a specified protocol.
Default: "illumina"

| Protocol      | Library Prep Kit                        | Trimming Parameter                   | 3' Adapter Sequence   |
| :------------ | :-------------------------------------- | :----------------------------------- | :-------------------  |
| zymo          | Zymo miRNA kit                          | clip_R1 = 1; three_prime_clip_R1 = 0 | TGGAATTCTCGGGTGCCAAGG |
| illumina      | Illumina TruSeq Small RNA               | clip_R1 = 0; three_prime_clip_R1 = 0 | TGGAATTCTCGGGTGCCAAGG |
| nextflex      | BIOO SCIENTIFIC  NEXTFLEX Small RNA-Seq | clip_R1 = 4; three_prime_clip_R1 = 4 | TGGAATTCTCGGGTGCCAAGG |
| qiaseq        | QIAGEN QIAseq miRNA                     | clip_R1 = 0; three_prime_clip_R1 = 0 | AACTGTAGGCACCATCAAT   |
| cats          | Diagenode CATS Small RNA-seq            | clip_R1 = 3; three_prime_clip_R1 = 0 | GATCGGAAGAGCACACGTCTG |

> Only 'zymo' and 'illumina' protocols have been tested, other protocols are inherited from nfcore/smrnaseq.

## Reference genomes
The pipeline config files come bundled with paths to the iGenomes reference index files.

### `--genome`
The reference genome to use of the analysis, needs to be one of the genome specified in the [config file](../conf/igenomes.config). This option is required. 

> Only 'Homo_sapiens(GRCh38)' and 'Rattus_norvegicus(Rnor_6.0)' have been tested for now.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the resources already set up. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      hairpin = '<path to the hairpin miRNA fasta file>'
      mirna_gtf = '<path to the miRNA gff file>' // mirbase GFF file
      mirtrace_species = "sps" // species according mirbase
      bowtie = '<path to the genome fasta file>' // Optional
      bed12 = '<path to the gene model BED12 file>' // Optional
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```
### `--mirna_gtf`
Use this to specify a miRNA GTF/GFF annotation file to use in the pipeline, instead of the one that comes with the selected genome.

### `--save_reference`
Supply this parameter to save any generated reference genome files to your results folder. These can then be used for future pipeline runs, reducing processing times.

## Trimming
### `--min_length`
The minimum length required for a read to be kept after trimming. Default: 18

### `--max_length`
The maximum length required for a read to be kept after trimming. Default: None

### `--save_trimmed`
Save trimmed sequences for future use. Default: false

### `--trim_nextseq`
Sets a quality score cutoff that ignores the quality scores of G bases. Default: 0 

## Skipping QC steps
### `--skip_fastqc`
Skip FastQC

### `--skip_multiqc`
Skip MultiQC

## Comparisons
### `--isomirs_fdr`
Use this to specify the false discovery rate (FDR) cutoff used in isomiRs. Default: 0.05

### `--isomirs_lfc`
Use this to specify the Log2 fold change cutoff in isomiRs. Please note that this is not a simple filter, this cutoff changes the null hypothesis in the statistical tests, therefore affects the p-values and FDRs. This is passed to DESeq2(isomiRs use DESeq2 in the backend) as a `Log2FoldChange` parameter, therefore, please use Log2 values of intended fold change cutoff. (Default is 0.585, which is Log2(1.5))

### `--comparisons`
By default, the pipeline will compare all pairwise combinations of sample groups. If this is not desirable, use this option to specify the path to a CSV file that describes which sample groups you want to compare. The CSV file should have the following format:
```
group_1,group_2
Experiment1,Control
Experiment2,Control
```
The header must be the same as shown. All group labels should also appear in the design CSV file "group" column.

## Job resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original). If it still fails after two times then the pipeline is stopped.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `awsbatch` profile and then specify following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch. This is required for `awsbatch` profile.

### `--awsregion`
The AWS region to run your job in. Default is set to `us-east-1`.

Please make sure to also set the `-work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters
### `-work-dir`
The working directory where intermediate files will be saved. This must be a S3 stroage path if using `awsbatch` profile.

### `--outdir`
The output directory where the results will be saved. This must be a S3 stroage path if using `awsbatch` profile.

### `--name`
Name used in report title and file name.

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

