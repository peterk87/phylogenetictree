#!/usr/bin/env python3

import click
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter



def count_ambig_nt(seq: str) -> int:
    counter = Counter(seq.lower())
    return sum(v for k,v in counter.items() if k not in {'a','g','c','t','-'})


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-lmin", "--lmin", help="ignore sequences w/length < lmin", required= False, type = int, default = 20000)
@click.option("-lmax", "--lmax", help="ignore sequences w/length >= lmax", required= False, type = int, default = 32000)
@click.option("-x", "--xambig", help="ignore sequences with >=  xambig ambiguous residues",required= False, type = int, default = 200)
@click.option("-i", "--gisaid_sequences", type=click.Path(exists=True), required= True)    
@click.option("-m", "--gisiad_metadata", type=click.Path(exists=True), required= True)
@click.option("-s", "--sample_lineage", type =str, required= True)
@click.option("-c", "--country", help="Country", required= False, type = str, default ='Canada')
@click.option("-r", "--region", help="Region", required= False, type= str, default = 'North America')
@click.option("-of", "--fasta_output", type=click.Path(exists=False), required= False)
@click.option("-om1", "--metadata_output1", help="New format of metadata", type=click.Path(exists=False), required= False)
@click.option("-om2", "--metadata_output2", help="Old format of metadata",type=click.Path(exists=False), required= False)

def main(lmin, lmax, xambig, gisaid_sequences, gisiad_metadata, sample_lineage, country, region, fasta_output, metadata_output1, metadata_output2):
     
    #Old format of GISAID Metadata, it is necessary to run further analysis nextstrain locally 

    column_names = ["strain", "virus", "gisaid_epi_isl", "genbank_accession", "date", "region", "country", "division", "location", "region_exposure", "country_exposure", "division_exposure", "segment", "length", "host", "age", "sex", "pangolin_lineage", "originating_lab", "submitting_lab", "authors", "url", "title", "date_submitted"]
    
    df_metadata = pd.DataFrame(columns = column_names)
    
    df = pd.read_table(gisiad_metadata)

    if (country != '' and region != ''):
        df_subset = df.loc[((df['Pango lineage'] == sample_lineage) & df['Location'].str.contains(region) & df['Location'].str.contains(country)), :]
    elif (country == '' and region != ''):
        df_subset = df.loc[((df['Pango lineage'] == sample_lineage) & df['Location'].str.contains(region)), :]
    elif (country != '' and region == ''):
        df_subset = df.loc[((df['Pango lineage'] == sample_lineage) & df['Location'].str.contains(country)), :]
    elif (country == '' and region == ''):
        df_subset = df.loc[((df['Pango lineage'] == sample_lineage)), :]

    strains_of_interest = set(df_subset['Virus name'])

    if len(strains_of_interest) > 0:
        with open(gisaid_sequences) as fin, open(fasta_output, 'w') as fout:
            for strains, seq in SimpleFastaParser(fin):

                if '|' in strains:
                    strains = strains.split('|')[0]

                if strains not in strains_of_interest:
                    continue
                
                if lmin < len(seq) <= lmax and count_ambig_nt(seq) < xambig:
                    if " " in strains:
                        strains.replace(" ","")
                    fout.write(f'>{strains}\n{seq}\n')
                    ########  Write  metadata to format as column_names #####
                    row = df_subset.loc[df['Virus name'] == strains]
                    strain = row['Virus name'].values[0]
                    virus ='ncov'
                    gisaid_epi_isl =row['Accession ID'].values[0]
                    genbank_accession ='?'
                    date = row['Collection date'].values[0]
                    found_location = row['Location'].values[0].split(' / ')
                    if (len(found_location) == 2):
                        region = found_location[0]
                        country= found_location[1]
                        division=''
                        location=''
                    elif (len(found_location) == 3):
                        region = found_location[0]
                        country= found_location[1]
                        division=found_location[2]
                        location=''
                    elif (len(found_location) == 4):
                        region = found_location[0]
                        country= found_location[1]
                        division=found_location[2]
                        location=found_location[3]
                    region_exposure = region
                    country_exposure = country
                    division_exposure = division
                    segment ='genome'
                    length = row['Sequence length'].values[0]
                    host = row['Host'].values[0]
                    age =row['Patient age'].values[0]
                    sex = row['Gender'].values[0]
                    pangolin_lineage = row['Pango lineage'].values[0]  
                    originating_lab=''
                    submitting_lab=''
                    authors=''
                    url=''
                    title=''
                    date_submitted = row['Submission date'].values[0]   
                    row_to_append = [strain, virus, gisaid_epi_isl, genbank_accession, date, region, country, division, location, region_exposure, country_exposure, division_exposure, segment, length, host, age, sex, pangolin_lineage, originating_lab, submitting_lab, authors, url, title, date_submitted]   
                    a_series = pd.Series(row_to_append, index = df_metadata.columns)
                    df_metadata = df_metadata.append(a_series, ignore_index=True)          
                else:
                    df_subset.drop(df_subset.loc[df['Virus name'] == strains].index, inplace = True)    
        # Keep the current format of metadata
        df_subset.to_csv(metadata_output1, index=False) 
        # Save to old format of metadata
        df_metadata.to_csv(metadata_output2, index=False)


if __name__ == '__main__':
	main()