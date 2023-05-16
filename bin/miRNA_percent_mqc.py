#!/usr/bin/env python

###
# This program calculates miRNA% from mirtrace results and output to the general stats table in MultiQC report
###

import argparse
import logging
import json
import csv
from collections import defaultdict

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def calc_miRNA_perc(results):

    # mirtrace results file is a TSV file
    with open(results, 'r') as fh:
        data = csv.DictReader(fh, delimiter="\t")
        # Get sample names
        rowname = data.fieldnames[0]
        snames = data.fieldnames[1:]
        if rowname != 'RNA_TYPE':
            logger.error("RNA_TYPE expected. Check if you have the right input file.")
        miRNA = dict()
        total = defaultdict(int)
        for row in data:
            # Find miRNA row
            if row[rowname] == 'miRNA':
                miRNA = dict(row)
                miRNA.pop(rowname)
                miRNA = { k:int(v) for k,v in miRNA.items() }
            # Calculate total
            for s in snames:
                total[s] += int(row[s])
        # Calculate percentage
        if len(miRNA):
            miRNA = { k:{'perc_miRNA':v/total[k]*100} for k,v in miRNA.items() }
        else:
            logger.error("No miRNA counts detected. Check your input.")
        
    # Base info for the general stats table
    gs_dict = {
        'id': 'perc_miRNA',
        'plot_type': 'generalstats',
        'pconfig': {
            'perc_miRNA': {
                'title': '% miRNA',
                'description': 'Percentage of miRNA reads in reads passing filter',
                'namespace': 'miRTrace',
                'scale': 'RdYlGn',
                'min': 0,
                'max': 100,
                'suffix': '%',
                'format': '{:,.1f}' 
            }
        }
    }

    # Add data to the table dict
    gs_dict['data'] = miRNA
    # Write the output to files
    with open('perc_miRNA_mirtrace_mqc.json', 'w') as ofh:
        json.dump(gs_dict, ofh, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse mirtrace results to calculate miRNA%""")
    parser.add_argument("mirtrace_rnatype_results", type=str, help="mirtrace rnatype results file")
    args = parser.parse_args()
    calc_miRNA_perc(args.mirtrace_rnatype_results)
    