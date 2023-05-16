#!/usr/bin/env python

###
# This program calculates the numbers of hairpin miRNA species detected from mirtrace results and output to the general stats table in MultiQC report
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

def calc_miRNA_species(results):

    # mirtrace results file is a TSV file
    with open(results, 'r') as fh:
        data = csv.DictReader(fh, delimiter="\t")
        # Get sample names
        rowname = data.fieldnames[0]
        snames = data.fieldnames[1:]
        if rowname != 'DISTINCT_MIRNA_HAIRPINS_ACCUMULATED_COUNT':
            logger.error("miRNA complexity results expected. Check if you have the right input file.")
        miRNA = defaultdict(int)
        for row in data:
            # If column of a sample is not empty, record the number of miRNA species, which is in first column
            mirna_species = int(row[rowname])
            for s in snames:
                if row[s] is not None and len(row[s]) > 0:
                    if mirna_species > miRNA[s]:
                        miRNA[s] = mirna_species
    
    # Covert the dict to correct format
    miRNA = { k:{'no_miRNA':v} for k,v in miRNA.items() }
        
    # Base info for the general stats table
    gs_dict = {
        'id': 'no_miRNA',
        'plot_type': 'generalstats',
        'pconfig': {
            'no_miRNA': {
                'title': 'No. miRNA',
                'description': 'Number of different miRNA hairpins detected in each sample',
                'namespace': 'miRTrace',
                'min': 0,
                'format': '{:,.0f}' 
            }
        }
    }

    # Add data to the table dict
    gs_dict['data'] = miRNA
    # Write the output to files
    with open('numbers_of_miRNA_mirtrace_mqc.json', 'w') as ofh:
        json.dump(gs_dict, ofh, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse mirtrace results to calculate no. miRNA detected""")
    parser.add_argument("mirtrace_mirna_complexity_results", type=str, help="mirtrace mirna complexity results file")
    args = parser.parse_args()
    calc_miRNA_species(args.mirtrace_mirna_complexity_results)
    