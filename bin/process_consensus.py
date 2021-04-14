#!/usr/bin/env python3
import collections
import sys
import re
import argparse
import logging
import os
from Bio import SeqIO
from os import listdir
from os.path import isfile, join, basename, abspath, isdir


def find_matched_files(input_paths,pattern,file_format='fasta',path_sep=',',folder_ignored=None):

    """Find files under given directories whose names contain given pattern. If
    input_paths inlcudes paths directing directly to files, then these files are added
    into output without the necessity of matching given pattern.
    """
    logger=logging.getLogger(__name__)       
    p=re.escape(pattern)
    if p!=pattern:
        logger.debug("pattern inlcudes regular expression metacharacters, escaped pattern is {}".format(p))
    matched_files=[]
    paths=list(map(abspath,input_paths.split(path_sep)))
    if folder_ignored is None: 
        folder_ignored=[]
    ol=len(paths)
    while(len(paths)):
        path=paths.pop(0)
        if isfile(path) and (ol>0 or re.search('^.*'+p+'.*'+file_format+'$',basename(path))):
            matched_files.append(path)            
        if isdir(path) and (basename(path) not in folder_ignored):
            paths.extend(map(lambda x: join(path,x),listdir(path)))  
        ol-=1
    logger.info('Files match with the identifier "{}": {}'.format(pattern,'\n'.join(matched_files)))
    return matched_files

def merge_consensus(args):
    
    matched_files = find_matched_files(args.consensus_path, args.identifier, file_format='', folder_ignored=['files'])
    matched_files.append(args.reference_fasta)

    completeness_threshold = float(args.completeness)
    with open (args.merge_consensus_sequence, 'w') as out_file:
        for files in matched_files:
            with open(files, 'r') as fasta_file:
                
                seq_records = SeqIO.parse(fasta_file, 'fasta')

                for record in seq_records:

                    #if (" " in record.id):
                        #record.id = record.id.replace(" ","_")

                    # require >75% completeness (non-N)
                    base_counter = collections.Counter()
                    for b in record.seq:
                        base_counter.update(b.upper())
                    total = 0
                    for base, count in base_counter.items():
                        total += count
                    completeness = 0
                    if total > 0:
                        completeness = 1 - (float(base_counter['N']) / total)
                    if completeness >= completeness_threshold:
                        SeqIO.write(record, out_file, 'fasta')
    
def main():
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-id","--identifier",help="Identifier for consensus sequences in fasta format",default='')
    parser.add_argument("-input","--consensus_path",help="Directory path of consensus sequencess",default='')
    parser.add_argument("-reference","--reference_fasta",help="Reference Sequence",default='')
    parser.add_argument("-output","--merge_consensus_sequence",help="Merge file of consensus sequence",default='')
    parser.add_argument('-c', '--completeness', default=0.75, help='the minimum completeness threshold (default: 0.75)')

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    merge_consensus(args)

if __name__ == "__main__":
    main()