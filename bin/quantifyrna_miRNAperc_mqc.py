#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import re
import json

def calc_perc_miRNA(RNAcounts):
    
    #Calculate percentage of miRNA as summed by quantify_rnacounts
    RNAcounts = pd.read_csv(RNAcounts, index_col = 0)
    perc_miRNA = round((RNAcounts["miRNA"]/RNAcounts.sum(axis=1))*100, 1)
    perc_miRNA = pd.DataFrame(perc_miRNA)
    perc_miRNA.columns = ["perc_miRNA_quantifyrna"]
    perc_miRNA = perc_miRNA.to_json(orient = "index")
    parsed = json.loads(perc_miRNA)

    gs_percmirna = {
            'id' : 'perc_miRNA_quantifyrna',
            'plot_type' : 'generalstats',
            'pconfig' : {
                'perc_miRNA_quantifyrna': {
                    'title': '% miRNA',
                    'description': 'Percentage of miRNA reads in reads passing filter',
                    'namespace' : 'Pipeline',
                    'scale' : 'RdYlGn',
                    'min' : 0,
                    'max' : 100,
                    'suffix' : '%',
                    'format' : '{:..1f}'
                    }
                }
            } 
    gs_percmirna['data'] = parsed

    with open('perc_miRNA_quantifyrna_mqc.json', 'w') as ofh:
        json.dump(gs_percmirna, ofh, indent = 4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Set up miRNA percentage to show in MultiQC General Stats Table""")
    parser.add_argument("-d", "--RNAtypecountdata", dest="RNAcounts", type=str, help="RNA Type quantifications compiled by quantify_rnacounts.py from rsem and tximport results")
    args = parser.parse_args()
    calc_perc_miRNA(args.RNAcounts)
