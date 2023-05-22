# Small RNA-Seq Nextflow Pipeline

### Introduction
This is a bioinformatics analysis pipeline used for small RNA sequencing data with a focus on miRNA. 

This pipeline was adapted from [nfcore/smrnaseq](https://github.com/nf-core/smrnaseq) version `1.0.0`. However, numerous changes have been made to the pipeline. They include, but are not limite to:
* Added differential miRNA expression analysis
* The library composition analysis now includes other small RNA types in addition to rRNA, tRNA, miRNAs
* Added differential expression analysis for genes of other small RNA types in addition to miRNAs
* Numerous report improvements including custom MultiQC modules

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

### Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Collapse reads ([`seqcsluter`](https://seqcluster.readthedocs.io/mirna_annotation.html#processing-of-reads))
4. Alignment against miRBase hairpin with reads from step 3 ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
5. Alignment with reads from step 2 against mature tRNA (GtRNAdb), mitochondrial tRNA (Ensembl), rRNA (Ensembl and UCSC Repeatmasker), lncRNA, scaRNA, snoRNA, snRNA, miscRNA (Ensembl), and miRNA hairpin (miRBase) ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml)) 
6. Read count quantification using RSEM with files from step 5 ([`RSEM`](http://deweylab.github.io/RSEM/README.html))
7. RSEM results summary with tximport ([`tximport`](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)) and RNA type quantification
8. Differential expression analysis of non-miRNA RNA types calculated with DESeq2 from step 7 tximport summary ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
9. miRNA and isomiR annotation from step 4 ([`mirtop`](https://github.com/miRTop/mirtop))
10. Sample comparison and statistical analysis ([`isomiRs`](https://www.bioconductor.org/packages/release/bioc/html/isomiRs.html))
11. miRNA quality control using reads from step 2 ([`mirtrace`](https://github.com/friedlanderlab/mirtrace))
12. Present stats, plots, results ([`MultiQC`](http://multiqc.info/))

### How to run the pipeline
It is preferable that you run the pipeline using [Aladdin Bioinformatics Platform](https://www.aladdin101.org/). Running this on your own computer or cloud is also possible. However, you would need to modify the genomes config file [conf/igenomes.config] to provide your own genome files.

#### Running on AWS Batch
``` bash
nextflow run Zymo-Research/aladdin-smrnaseq \
    --design "<design CSV file>" \
    --genome 'Homo_sapiens(GRCh38)' \
    --protocol zymo \
    -profile awsbatch \
    -work-dir "<work dir on S3>" \
    --awsqueue "<SQS ARN>" \
    --outdir "<output dir on S3>" \
    --name "my_analysis" \
    -r "0.0.1" \
```
1. The options `-design` and `--genome` are required. Remember to enclose `--genome` parameter value in () because it contains special characters.
2. The options `-profile`, `-work-dir`, `--outdir`, and `--awsqueue` are required only when running pipelines on AWS Batch, but are highly recommended.
3. The option `-r` helps pin workflow to a specific release on GitHub.
5. The option `--name` will add a custom title to the report.
6. The design CSV file must have the following format.
```
group,sample,read_1,read_2
Control,Sample1,s3://mybucket/this_is_s1_R1.fastq.gz,
Control,Sample2,s3://mybucket/this_is_s2_R1.fastq.gz,
Experiment,Sample3,s3://mybucket/that_is_s3_R1.fastq.gz,
Experiment,Sample4,s3://mybucket/that_be_s4_R1.fastq.gz,
```
The header line must be present and cannot be changed. Sample labels and group names must contain only alphanumerical characters or underscores, must not start with "R1" or "R2", and must start with a letter. Sample labels and group names also cannot include phrases that will be automatically filtered by MultiQC. A list of terms unsuitable for sample and group labels in this pipeline can be viewed in the MultiQC source code [here](https://github.com/ewels/MultiQC/blob/b936a7a6d7050f3edc1ceefe8ae6ecd93865bf66/multiqc/utils/config_defaults.yaml#L150-L284). Only single-end sequencing is supported. Full S3 paths of R1 FASTQ files should be provided. R2 FASTQ files are ignored if provided. If you do not wish to carry out statistical comparisons of samples, simply leave the group column blank (but keep the comma).

#### Running with Docker
``` bash
nextflow run Zymo-Research/aladdin-smrnaseq \
    --design "<design CSV file>" \
    --genome 'Homo_sapiens(GRCh38)' \
    --protocol zymo \
    -profile docker \
    --outdir "<output dir on S3>" \
    --name "my_analysis" \
    -r "0.0.1" \
```

For more options to customize your run, please refer to [usage documentation](docs/usage.md).
