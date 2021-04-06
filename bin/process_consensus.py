#!/usr/bin/env python3

import collections
import pysam
import sys
import argparse


def merge_consensus(args):
    print (args.global_lineage_report)
            

    
def main():
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-g','--global_lineage_report', required=True, help='Lineage report of global samples from GISIAD')
    parser.add_argument('-c', '--samples_lineage_report', required=True, help='Lineage report of samples')
    parser.add_argument('-o', '--samples_lineage_report', required=True, help='Lineage report of samples')

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    merge_consensus(args)

if __name__ == "__main__":
    main()