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
        '4': {'start': 2e5, 'stop': 9e5, 'buffer': 2e4},
        'X': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        'Y': {'start': 1e6, 'stop': 2e6, 'buffer': 1e5}
    }

    tumour_svs = {
        1: [200, 'deletion', 'None', 0],
        2: [100, 'deletion', 'None', 0],
        3: [20, 'tandem duplication', 1, 0],
        4: [5, 'inversion', 'None', 0],
        5: [300, 'deletion', 'None', 0],
        6: [2, 'deletion', 'None', 0],
        7: [100, 'tandem duplication', 1, 0],
        8: [50, 'deletion', 'None', 0],
        9: [25, 'deletion', 'None', 0],
        10: [1, 'inversion', 'None', 0]
    }

    normal_svs = {
        3: '2L',
        6: '3R',
    }

    chroms = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']

    return chroms, sv_regions, tumour_svs, normal_svs

def parse_bed(options):

    chroms, sv_regions, tumour_svs, normal_svs = set_regions()

    th1, th2 = 'tumour_h1.bed', 'tumour_h2.bed'
    nh1, nh2 = 'normal_h1.bed', 'normal_h2.bed'


    with open(th1, 'w') as t1, open(th2, 'w') as t2, open(nh1, 'w') as n1, open(nh2, 'w') as n2:
        for chrom in chroms:
            region_start, window_end = int(sv_regions[chrom]['start']), int(sv_regions[chrom]['stop'])
            sv_start = region_start
            for sv in tumour_svs.keys():
                if sv > options.svs_per_chrom: continue
                sv_length = tumour_svs[sv][0] * 1000
                items = '\t'.join(map(str, tumour_svs[sv][1:]))
                sv_end = sv_start + sv_length
                if sv_end >= window_end:
                    print("The end of this SV %s is bigger than the region on chrom %s!" % (tumour_svs[sv], chrom))
                l = '\t'.join(map(str, [chrom, sv_start, sv_end, items]))

                write_haplotypes(l, sv, t1, t2)
                if sv in normal_svs and normal_svs[sv] == chrom:
                    write_haplotypes(l, sv, n1, n2)

                sv_start += sv_length + int(sv_regions[chrom]['buffer'])

    return True


def write_haplotypes(l, sv, h1, h2):
    if sv % 2 != 0:
        h1.write(l + '\n')
    else:
        h2.write(l + '\n')


def get_args():
    parser = OptionParser()

    parser.add_option("-d", "--tumour-depth", dest="depth", action="store", help="Average sequencing depth for tumour sample. [10]")
    parser.add_option("-n", "--svs-per-chrom", dest="svs_per_chrom", action="store", help="Number of SVs per chromosomes. [10]")

    parser.set_defaults(bedfile='svs.bed', depth=30, svs_per_chrom=10)

    options, args = parser.parse_args()

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