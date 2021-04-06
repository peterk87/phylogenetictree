#!/usr/bin/env python3

import os
import logging
import argparse
from Bio import SeqIO
import pandas as pd
from collections import Counter


def usage_description():

    des="Filter global sequences based on length and number of ambiguous positions, parent lineage, the script will output new fasta file\n"
    des="\n"+des+"\n"    
    return des

def init_config():

    parser = argparse.ArgumentParser(description=usage_description()) 

    parser.add_argument("--lmin",help="remove sequences w/length < lmin", required= True)
    parser.add_argument("--xambig",help="remove sequences with >=  xambig ambiguous residues",required= True)
    parser.add_argument("--fasta_input",help="Fasta sequences after cleaning", required= True)    
    parser.add_argument("--metadata",help="Metadata from GISIAD", required= True)
    parser.add_argument("--sample_lineage",help="Pangolin Lineage of sample", required= False)
    parser.add_argument("--country",help="Country", required= False)
    parser.add_argument("--region",help="Region", required= False)
    parser.add_argument("--fasta_output",help="Fasta sequences after filtering", required= True)
    parser.add_argument("-d","--debug",help="Debug Mode",action="store_true") 
    
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    return args

def count_ambig_nt(seq: str) -> int:
    counter = Counter(seq.lower())
    return sum(v for k,v in counter.items() if k not in {'a','g','c','t','-'})


def main(args):

    ########### Get metadata information #################################
    
    df_meta = pd.read_table(args.metadata)

    ########### Grab sequences that belong to parent lineages ############
    lineage_search =''
    if '.' in args.sample_lineage:
        if len(args.sample_lineage.split('.')) > 2:
            lineage_search = args.sample_lineage[:3]
        else:
            lineage_search = args.sample_lineage
    else:
        lineage_search = args.sample_lineage

    df_subset = df_meta.loc[((df_meta.country == args.country) & (df_meta.pangolin_lineage == args.sample_lineage)), :]

    print (df_subset)

if __name__=='__main__':
    main(init_config())