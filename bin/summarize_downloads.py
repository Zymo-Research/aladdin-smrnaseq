#!/usr/bin/env python
###
# This program creates a JSON file that lists the locations of files for download
# and how to display those files on aladdin platform.
###

import argparse
import logging
import os
import json
import csv

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def summarize_downloads(locations, design):

    """
    :param locations: a file containing locations of files on S3
    :param design: the design file containing group and sample labels
    """

    file_info = dict()

    # Read the design file to collect valid sample and group labels
    logger.info("Reding design file...")
    groups = set()
    samples = set()
    with open(design, 'r') as fh:
        data = csv.DictReader(fh)
        for row in data:
            samples.add(row['sample'])
            if len(row['group']):
                groups.add(row['group'])

    # Define what to do with each type of files
    categories = {
        '_multiqc_report.html': ('Report', 'report'),
        '_isomiRs_MA_plot.jpg': ('MA plot of DE miRNA', 'comparisons'),
        'isomiRs_sample_MDS_plot.jpg': ('miRNA MDS plot of samples', 'all_samples'),
        '_isomiRs_scatterplot.jpg': ('miRNA Scatterplot', 'comparisons'),
        'isomiRs_sample_similarity_matrix.jpg': ('miRNA sample similarity matrix', 'all_samples'),
        'isomiRs_top_miRNA_gene_heatmap.jpg': ('Heatmap of top DE miRNA', 'all_samples'),
        '_isomiRs_results.xlsx': ('miRNA DE Analysis Results', 'comparisons'),
        '_deseq2_MA_plot.jpg': ('MA plot of DE non-miRNA', 'comparisons'),
        'deseq2_sample_MDS_plot.jpg': ('non-miRNA MDS plot of samples', 'all_samples'),
        '_deseq2_scatterplot.jpg': ('non-miRNA Scatterplot', 'comparisons'),
        'deseq2_sample_similarity_matrix.jpg': ('non-miRNA sample similarity matrix', 'all_samples'),
        'deseq2_top_nonmiRNA_gene_heatmap.jpg' : ('Heatmap of top DE non-miRNA genes', 'all_samples'),
        '_deseq2_results.xlsx' : ('non-miRNA DEG analysis results', 'comparisons')
    }

    # Read the file locations
    with open(locations, 'r') as fh:
        for line in fh:
            info = dict()
            path = line.strip()
            info['path'] = path
            logger.info("Processing file {}".format(path))
            # Get file name
            fn = os.path.basename(path)
            # Check each file category
            for suffix, values in categories.items():
                if fn.endswith(suffix):
                    file_type, scope = values
                    info['file_type'] = file_type
                    info['scope'] = scope
                    if scope in ['samples', 'comparisons']:
                        sname = fn.replace(suffix, '')
                        if scope == 'samples':
                            # Check if the parsed sample name is in the original design
                            if sname in samples:
                                info['sample'] = sname
                            else:
                                logger.error("Parsed sample name from {} not found in the design file".format(fn))
                        else:
                            # Check if the parsed group names are in the original design
                            g1, g2 = sname.split('_vs_')
                            if g1 in groups and g2 in groups:
                                info['comparison'] = sname
                            else:
                                logger.error("Parsed group names from {} not found in the design file".format(fn))
                    file_info[fn] = info
                    break
            else:
                logger.info("File {} did not match any expected patterns".format(fn))
    
    # Output the dict to JSON
    with open('files_to_download.json', 'w') as fh:
        json.dump(file_info, fh, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Generate a json file for displaying outputs on aladdin platform""")
    parser.add_argument("locations", type=str, help="File with all the locations of files on S3")
    parser.add_argument("-d", "--design", dest="design", required=True, type=str, help="Design file for sanity check purposes")
    args = parser.parse_args()
    summarize_downloads(args.locations, args.design)  
