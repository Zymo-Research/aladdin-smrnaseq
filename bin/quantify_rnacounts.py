#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import pysam
import re
import json

def sumtximport(tximporttsv, idmaptorna, bamlist, rsemlogs):
    
    #import tximport summary
    tximporttsv = pd.read_csv(tximporttsv, sep="\t")
    #idmaptorna contains information assigning FASTA IDs to rna type
    idmaptorna = pd.read_csv(idmaptorna, sep = "\t")
    idmaptorna = idmaptorna.set_index('transcript_id')
    RNAcounts = pd.DataFrame()
    RNAtypes = ["miRNA", "tRNA", "rRNA", "lncRNA", "miscRNA", "scaRNA", "snoRNA", "snRNA" ]

    #label FASTA IDs in tximport tsv file with RNA types as recorded in the ID map
    #Sum total of each RNA type and concatenate to RNAcounts table
    for RNAtype in RNAtypes:
        RNAtypelist = idmaptorna[idmaptorna["transcript_biotype"] == RNAtype]
        RNAtypecounts = RNAtypelist.join(tximporttsv, how = "inner").drop(["transcript_biotype"], axis=1)
        RNAsum = RNAtypecounts.sum(axis=0)
        RNAcounts = pd.concat((RNAcounts, RNAsum), axis=1)

    RNAcounts.columns = RNAtypes
    #Change file names from RNAcounts csv to only sample names; sort rows based on sample names
    indexreplace = RNAcounts.index.str.replace(".isoforms.results", "", regex=False)
    RNAcounts.index = indexreplace
    RNAcounts = RNAcounts.sort_index()

    #Create dataframes which will count unmapped reads using BAM files, and number of reads dropped by rsem
    unmappeddf = pd.DataFrame()
    thrownreadsdf = pd.DataFrame()
    totalreadspattern = "Thread.*: N = (\d+)"
    
    for index in RNAcounts.index:
        #Extract number of reads thrown out by rsem
        totalreads = 0
        with open(index+"_rsem.log", "r") as rsem_log:
            for line in rsem_log:
                if re.search(totalreadspattern, line):
                    pattern_set = re.compile(totalreadspattern)
                    match = pattern_set.match(line)
                    totalreads += float(match.group(1))

        thrownreads = totalreads - RNAcounts.sum(axis=1)[index]
        thrownreads = pd.DataFrame(thrownreads, index=[index], columns = ["excluded_by_rsem"])
        thrownreadsdf = pd.concat([thrownreadsdf, thrownreads], axis = 0)
        #Extract number of reads unmapped by bowtie
        bamname = index+".bam"
        unmapped = pysam.view("-c", "-f", "4", bamname)
        unmapped = float(re.sub("\n", "", unmapped))
        unmapped = pd.DataFrame(unmapped, index=[index], columns = ["unmapped"]) 
        unmappeddf = pd.concat([unmappeddf, unmapped], axis = 0)
    
    RNAcounts = pd.concat([RNAcounts, thrownreadsdf,  unmappeddf], axis = 1)
    RNAcounts.to_csv("RNAcounts.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Quantify read numbers in bam files""")
    parser.add_argument("-t", "--tximportdata", dest="tximporttsv", type=str, help="quantifications compiled by tximport from rsem results")
    parser.add_argument("-f", "--bamfiles", dest="bamlist", type=str,  help="list of bam files prepared by bowtie")
    parser.add_argument("-i", "--idmaptornatype", dest="idmaptorna", type=str, help="file matching gtrnadb, mitornadb, ensembl, repeatmasker rRNA, miRBase codes to RNA types")
    parser.add_argument("-r", "--rsemlogs", dest="rsemlogs", type=str, help="log files produced by rsem")
    args = parser.parse_args()
    sumtximport(args.tximporttsv, args.idmaptorna, args.bamlist, args.rsemlogs)
