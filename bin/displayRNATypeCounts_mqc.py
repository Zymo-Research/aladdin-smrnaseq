#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import re
import json

def displayRNATypeCounts(RNAcounts, max_length):

    #Set json file for MultiQC RNA Type Counts Chart
    RNAcounts = pd.read_csv(RNAcounts, index_col = 0)
    RNAcounts_json = RNAcounts.to_json(orient = "index")
    RNAcounts_parsed = json.loads(RNAcounts_json)

    description = "The following bargraph provides an estimate of RNA types in the sample. Reads were mapped to RNA reference sequences with Bowtie, and subsequent transcript quantification was completed with RSEM. Current piRNA database entries may contain significant sequential overlap with miRNA entries and negatively affect miRNA results. To preserve miRNA results, piRNA were not included in this analysis. \"Unmapped\" reads may include RNA types that the pipeline does not test for, such as circRNA and piRNA. \"miscellaneous\" RNA, or miscRNA, includes unclassified RNA types, such as vault RNA and YRNA."
    maxlength_message = " The current Trim Galore max length setting filters out all reads longer than " + max_length +" basepairs. Longer RNA types may be undercounted."
    if(max_length != ""):
        description = description + maxlength_message

    RNAtypechart_multiqc = {
            'id' : 'RNAtypechart',
            'section_name' : 'Estimated RNA Type Counts',
            'description' : description,
            'plot_type' : 'bargraph',
            'pconfig' : {
                    'id': 'estimated_RNA_type_counts_plot',
                    'title': 'RNA Type Plot',
                    'ylab': 'Counts'
                    },
            'categories' : {
                    'miRNA': {
                        'name': 'miRNA',
                        'color': '#3A4BA2'
                        },
                    'tRNA': {
                        'name': 'tRNA',
                        'color': '#3EE4BF'
                        },
                    'rRNA': {
                        'name': 'rRNA',
                        'color': '#A7122A'
                        },
                    'lncRNA': {
                        'name': 'lncRNA',
                        'color': '#ECE343'
                        },
                    'miscRNA' : {
                        'name': 'miscRNA',
                        'color': '#0CAD0E'
                        },
                    'scaRNA': {
                        'name': 'scaRNA',
                        'color': '#E43DCA'
                        },
                    'snoRNA': {
                        'name': 'snoRNA',
                        'color': '#ECB289'
                        },
                    'snRNA': {
                        'name': 'snRNA',
                        'color': '#ADD7F7'
                        },
                    'excluded_by_rsem' : {
                        'name': 'excluded_by_rsem',
                        'color': '#123456'
                        },
                     'unmapped' : {
                        'name': 'unmapped',
                        'color': '#999999'
                        }
                    }
        }

    RNAtypechart_multiqc['data'] = RNAcounts_parsed
    with open('RNAtypechart_mqc.json', 'w') as ofh:
         json.dump(RNAtypechart_multiqc, ofh, indent = 4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Set up quantified RNA type counts to show in MultiQC Chart""")
    parser.add_argument("-d", "--RNAtypecountdata", dest="RNAcounts", type=str, help="RNA Type quantifications compiled by quantify_rnacounts.py from rsem and tximport results")
    parser.add_argument("-m", "--max_length", dest = "max_length", type=str, help="User assigned maximum read length of the current pipeline run")
    args = parser.parse_args()
    displayRNATypeCounts(args.RNAcounts, args.max_length)
