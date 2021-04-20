#!/usr/bin/env python3

import click
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter



def count_ambig_nt(seq: str) -> int:
    counter = Counter(seq.lower())
    return sum(v for k,v in counter.items() if k not in {'a','g','c','t','-'})


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-lmin", "--lmin", help="remove sequences w/length < lmin", required= False, type = int, default = 26000)
@click.option("-lmax", "--lmax", help="remove sequences w/length >= lmax", required= False, type = int, default = 32000)
@click.option("-x", "--xambig", help="remove sequences with >=  xambig ambiguous residues",required= False, type = int, default = 200)
@click.option("-i", "--gisaid_sequences", type=click.Path(exists=True), required= True)    
@click.option("-m", "--gisiad_metadata", type=click.Path(exists=True), required= True)
@click.option("-s", "--sample_lineage", required= False)
@click.option("-c", "--country", help="Country", required= False, default ='Canada')
@click.option("-r", "--region", help="Region", required= False)
@click.option("-of", "--fasta_output", type=click.Path(exists=False), required= False)
@click.option("-om", "--metadata_output", type=click.Path(exists=False), required= False)

def main(lmin, lmax, xambig, gisaid_sequences, gisiad_metadata, sample_lineage, country, region, fasta_output, metadata_output):

    pangolin_lineages = ['B.1.1.306']

    df = pd.read_table(gisiad_metadata)
    
    df_subset = df.loc[((df.pangolin_lineage.isin(pangolin_lineages)) & (df.country == country)), :]

    df_subset.to_csv(metadata_output, index=False)  

    strains_of_interest = set(df_subset.strain)

    with open(gisaid_sequences) as fin, open(fasta_output, 'w') as fout:
        for strains, seq in SimpleFastaParser(fin):
            if strains not in strains_of_interest:
                continue
            if lmin < len(seq) <= lmax and count_ambig_nt(seq) < xambig:
                fout.write(f'>{strains}\n{seq}\n')


if __name__ == '__main__':
	main()