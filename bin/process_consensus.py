#!/usr/bin/env python3
import collections
import sys
import re
import argparse
import logging
import os
from Bio import SeqIO
from os import listdir
#from os.path import isfile, join, basename, abspath, isdir
from pathlib import Path, PurePath

logger=logging.getLogger(__name__) 

def find_matched_files(input_paths,pattern,file_format='fasta',path_sep=',',folder_ignored=None):

    """Find files under given directories whose names contain given pattern. If
    input_paths inlcudes paths directing directly to files, then these files are added
    into output without the necessity of matching given pattern.
    """
    p=re.escape(pattern)
    if p!=pattern:
        logger.debug("pattern inlcudes regular expression metacharacters, escaped pattern is {}".format(p))
    matched_files=[]
    paths = [Path(pathx).resolve() for pathx in input_paths.split(path_sep)]
    if folder_ignored is None: 
        folder_ignored=[]
    ol=len(paths)
    while(len(paths)):
        path=paths.pop(0)
        if Path(path).is_file() and (ol>0 or re.search('^.*'+p+'.*'+file_format+'$', PurePath(path).name)):
            matched_files.append(path)            
        if Path(path).is_dir() and (PurePath(path).name not in folder_ignored):
            paths.extend(map(lambda x: PurePath(path).joinpath(x), Path(path).iterdir()))  
        ol-=1
    logger.info('Files match with the identifier "{}": {}'.format(pattern, matched_files))
    return matched_files

def merge_consensus(args):
    
    matched_files = find_matched_files(args.consensus_path, args.identifier, file_format='', folder_ignored=['files'])
    completeness_threshold = float(args.completeness)
    with open (args.output, 'w') as out_file:
        for files in matched_files:
            with open(files, 'r') as fasta_file:
                
                seq_records = SeqIO.parse(fasta_file, 'fasta')

                for record in seq_records:
                    # require >75% completeness (non-N)
                    if args.filter_incomplete_sequences:
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
                    else:
                        SeqIO.write(record, out_file, 'fasta')
    
def main():
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-id","--identifier",help="Identifier for consensus sequences in fasta format",default='')
    parser.add_argument("-i","--consensus_path",help="Directory path of consensus sequencess",default='')
    parser.add_argument("-o","--output",help="Merge file of consensus sequence",default='')
    parser.add_argument("-f","--filter_incomplete_sequences",help="Filter incomplete sequences",default=False)
    parser.add_argument('-c', '--completeness', help='the minimum completeness threshold (default: 0.75)', default=0.75)

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    merge_consensus(args)

if __name__ == "__main__":
    main()