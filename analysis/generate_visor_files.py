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

def set_regions(options):

    sv_regions = {
        '2L': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '2R': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '3L': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '3R': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        '4': {'start': 2e5, 'stop': 10e5, 'buffer': 2e4},
        'X': {'start': 5e6, 'stop': 15e6, 'buffer': 1e6},
        'Y': {'start': 7e5, 'stop': 2.5e6, 'buffer': 1e5}
    }

    genome_regions = {
        '2L': {'start': 1e6, 'stop': 22e6},
        '2R': {'start': 1e6, 'stop': 25e6},
        '3L': {'start': 1e6, 'stop': 28e6},
        '3R': {'start': 1e6, 'stop': 30e6},
        '4': {'start': 1e5, 'stop': 12.5e5},
        'X': {'start': 1e6, 'stop': 22e6},
        'Y': {'start': 2e5, 'stop': 3.5e6}
    }

    purities = [100.0, 80.0, 60.0, 40.0, 20.0]
    coverage_flux = 80.0

    write_regions(genome_regions, purities, coverage_flux, options)

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


def write_regions(genome_regions, purities, coverage_flux, options):
    for p in purities:
        tumour_regions = os.path.join(options.out_dir, 'tumour_regions_' + str(int(coverage_flux)) + '_' + str(int(p)) + '.bed')
        normal_regions = os.path.join(options.out_dir, 'normal_regions_' + str(int(coverage_flux)) + '_' + str(int(p)) + '.bed')
        with open(tumour_regions, 'w') as tr, open(normal_regions, 'w') as nr:
            for c in genome_regions.keys():
                start = int(genome_regions[c]['start'])
                end = int(genome_regions[c]['stop'])
                t = '\t'.join(map(str, [c, start, end, coverage_flux, p]))
                n = '\t'.join(map(str, [c, start, end, coverage_flux, 100.0]))
                tr.write(t + '\n')
                nr.write(n + '\n')

    return True


def parse_bed(options):

    chroms, sv_regions, tumour_svs, normal_svs = set_regions(options)

    th1 = os.path.join(options.out_dir, 'tumour_h1.bed')
    th2 = os.path.join(options.out_dir, 'tumour_h2.bed')
    nh1 = os.path.join(options.out_dir, 'normal_h1.bed')
    nh2 = os.path.join(options.out_dir, 'normal_h2.bed')

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
                    print("The end [%s] of this SV %s occurs after the end [%s] of chrom %s!" % (sv_end, tumour_svs[sv], window_end, chrom))
                    continue
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
    parser.add_option("-o", "--out-dir", dest="out_dir", action="store", help="Directory to write bed files to [.]")

    parser.set_defaults(bedfile='svs.bed', depth=30, svs_per_chrom=10, out_dir='.')

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