#!/usr/bin/env python

import argparse
import re
import os

def check_design(DesignFileIn, DesignFileOut, Comparison):

    HEADER = ['group','sample','read_1','read_2']
    legal_pattern = r"^[a-zA-Z][a-zA-Z0-9_]*$"
    illegal_patterns = ["_tophat", "ReadsPerGene", "_star_aligned", "_fastqc", "_counts", "Aligned", "_slamdunk", "_bismark", "_SummaryStatistics", \
                        "_duprate", "_vep", "ccs", "_NanoStats", "_R1", "_R2", "_trimmed", "_val", "_mqc", "short_summary_", "_summary", "_matrix", \
                        r"_$", r"^R1", r"^R2", "_isomiRs_results", "_deseq2_results", "_dinucleotide_frequency", "_rseqc"]
    
    with open(DesignFileIn, "r") as fin:
        # Check the header first
        header = fin.readline().strip().split(',')
        assert header == HEADER, "Header of design file Incorrect! Should be {}".format(','.join(HEADER))
        
        labels = []
        groups = []
        # Check the rest
        for line in fin:
            cols = line.strip().split(',')
            # Check the number of columns in each line
            assert len(cols) == len(HEADER), "Number of columns incorrect in line '{}'!".format(line.strip())
            # Check the sample label
            assert re.match(legal_pattern, cols[1]), "Sample label {} contains illegal characters or does not start with letters!".format(cols[1])
            for phrase in illegal_patterns:
                assert re.search(phrase, cols[1]) == None, "Sample label {} contains file phrase(s) that will be automatically filtered out by MultiQC in the final report. Please choose a different label.".format(cols[1])
            assert cols[1] not in labels, "Duplicate sample label {}".format(cols[1])
            labels.append(cols[1])
            # Check the fastq file locations
            assert cols[2], "R1 path could not be missing!"
            assert cols[2].endswith('fastq.gz') or cols[2].endswith('fq.gz'), "R1 path must end with fastq.gz or fq.gz!"
            if cols[3]:
                print("R2 ignored in smRNAseq.")
            # Check the group label
            if cols[0]:
                assert re.match(legal_pattern, cols[0]), "Group label {} contains illegal characters or does not start with letters!".format(cols[0])
                for phrase in illegal_patterns:
                    assert re.search(phrase, cols[0]) == None, "Group label {} contains file phrase(s) that may be automatically filtered out by MultiQC in the final report. Please choose a different label.".format(cols[0])
                groups.append(cols[0])
            assert len(groups)==0 or len(groups) == len(labels), "Group label(s) missing in some but not all samples!"

    if Comparison is not None:
        HEADER = ['group_1', 'group_2']
        with open(Comparison, "r") as fin:
            # Check the header first
            header = fin.readline().strip().split(',')
            assert header == HEADER, "Header of group comparison file Incorrect! Should be {}".format(','.join(HEADER))
            # Check the groups
            for line in fin:
                cols = line.strip().split(',')
                # Check the number of columns in each line
                assert len(cols) == len(HEADER), "Number of columns incorrect in line '{}'!".format(line.strip())
                # Check the group labels
                assert cols[0] in groups, "Group label '{}' not in design file".format(cols[0])
                assert cols[1] in groups, "Group label '{}' not in design file".format(cols[1])
    
    os.rename(DesignFileIn, DesignFileOut)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Sanity check the smRNAseq design CSV file and group comparison file""")
    parser.add_argument("DesignFileIn", type=str, help="Input design CSV file")
    parser.add_argument("-c", "--comparison", dest="Comparison", help="Group comparison CSV file")
    args = parser.parse_args()
    check_design(args.DesignFileIn, "checked_"+args.DesignFileIn, args.Comparison)
