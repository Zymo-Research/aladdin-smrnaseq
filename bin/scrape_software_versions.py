#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'smrnaseq pipeline': ['v_pipeline.txt', r"(\S+)"],
    'R': ['v_R.txt', r"R version (\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'Bowtie': ['v_bowtie.txt', r"version (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Seqcluster': ['v_seqcluster', r"seqcluster\s+(\S+)"],
    'FASTX': ['v_fastx.txt', r"Toolkit (\S+)"],
    'mirtop' : ['v_mirtop.txt', r"mirtop\s+(\S+)"],
    'miRTrace': ['v_mirtrace.txt', r"(\S+)"],
    'rsem': ['v_rsem.txt', r"RSEM v(\S+)"],
    'tximport': ['v_tximport.txt', r"(\S+)"],
    'isomiRs': ['v_isomiRs.txt', r"(\S+)"],
    'DESeq2' : ['v_DESeq2.txt', r"(\S+)"],
}
results = OrderedDict()
results['smrnaseq pipeline'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['Seqcluster'] = '<span style="color:#999999;\">N/A</span>'
results['FASTX'] = '<span style="color:#999999;\">N/A</span>'
results['miRTrace'] = '<span style="color:#999999;\">N/A</span>'
results['rsem'] = '<span style="color:#999999;\">N/A</span>'
results['tximport'] = '<span style="color:#999999;\">N/A</span>'
results['mirtop'] = False
results['isomiRs'] = False
results['DESeq2'] = False


# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'Software Versions'
plot_type: 'html'
description: 'Software versions are collected at run time from the software output. This pipeline is adapted from <a href="https://github.com/nf-core/smrnaseq" target="_blank">nf-core smRNAseq pipeline</a>.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    # Only display versions of softwares that were actually used.
    if v:
        print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        if v:
            f.write("{}\t{}\n".format(k,v))
