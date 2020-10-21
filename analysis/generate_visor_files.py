from __future__ import division
from __future__ import print_function

import fnmatch
import json
import ntpath
import os, re
import sys
from collections import defaultdict
from optparse import OptionParser

import pandas as pd

def set_regions():

    sv_regions = {
        '2L': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '2R': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '3L': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '3R': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '4': {'start': 4e5, 'stop': 9e5, 'buffer': 2e4},
        'X': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        'Y': {'start': 1e6, 'stop': 2e6, 'buffer': 1e5}
    }

    svs_to_generate = {
        1: [50, 'deletion', 'None', 0],
        2: [100, 'deletion', 'None', 0],
        3: [20, 'tandem duplication', 1, 0],
        4: [5, 'inversion', 'None', 0],
        5: [10, 'deletion', 'None', 0],
        6: [50, 'deletion', 'None', 0],
        7: [100, 'tandem duplication', 1, 0],
        8: [20, 'deletion', 'None', 0],
        9: [5, 'deletion', 'None', 0],
        10: [1, 'inversion', 'None', 0]
    }

    chroms = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']

    return chroms, sv_regions, svs_to_generate

def parse_bed(options):

    svs_per_chrom = 5
    chroms, sv_regions, svs_to_generate = set_regions()


    with open(options.bedfile, 'w') as bed_out:
        for chrom in chroms:
            start, end = int(sv_regions[chrom]['start']), int(sv_regions[chrom]['stop'])
            walker = start
            length = end - start
            for sv in svs_to_generate.keys():
                sv_length = svs_to_generate[sv][0] * 1000
                items = '\t'.join(map(str, svs_to_generate[sv][1:]))
                l = '\t'.join(map(str, [chrom, walker, walker + sv_length, items]))
                # print(l)
                bed_out.write(l + '\n')
                walker += sv_length + int(sv_regions[chrom]['buffer'])



    return True


def get_args():
    parser = OptionParser()

    parser.add_option("-b", "--bedfile", dest="bedfile", action="store", help=".bed file to write SVs to")
    # parser.add_option("-g", "--genome", dest="genome", action="store", help="Genome fasta file")
    # parser.add_option("-w", "--window", dest="window", action="store", type="int", help="Window to slide over sequence. [Default 5]")
    # parser.add_option("-l", "--locus", dest="locus", action="store", help="chr:start-end")
    # parser.add_option("--event", dest="event", action="store", help="Inspect one event")

    parser.set_defaults(bedfile='svs.bed')

    options, args = parser.parse_args()

    if not options.bedfile:
        parser.print_help()
        print()
        sys.exit("[!] Must provide a bed file. Exiting.")
    else:
        return options, args


def main():
    options, args = get_args()
    try:
        parse_bed(options)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return


if __name__ == "__main__":
    sys.exit(main())