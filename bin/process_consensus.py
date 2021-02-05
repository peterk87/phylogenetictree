#!/usr/bin/env python3

import collections
import pysam
import sys
import argparse


def merge_consensus(args):
    
    completeness_threshold = float(args.completeness)  

    for fn in args.consensusfiles:

        with pysam.FastxFile(fn) as fa:
            for record in fa:
                # shrink the long name ivar uses
                #record.name = record.name.replace(".primertrimmed.consensus_threshold_0.75_quality_20", "")

                # require >75% completeness (non-N)
                
                base_counter = collections.Counter()
                for b in record.sequence:
                    base_counter.update(b.upper())

                total = 0
                for base, count in base_counter.items():
                    total += count
                completeness = 0
                if total > 0:
                    completeness = 1 - (float(base_counter['N']) / total)
                if completeness >= completeness_threshold:
                    sys.stdout.write(">" + record.name + "\n")
                    sys.stdout.write(record.sequence + "\n")
            

    
def main():
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--consensusfiles', required=True)
    parser.add_argument('-c', '--completeness', default=0.75,
                    help='the minimum completeness threshold (default: 0.75)')

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    merge_consensus(args)

if __name__ == "__main__":
    main()