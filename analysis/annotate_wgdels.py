#!/usr/bin/env python
import os, re, sys
import ntpath
import pandas as pd
import csv

wgdels_in = 'data/wholegut_dels_merged.bed'
wgvars_in = 'data/all_WG_samples_merged_filt.txt'

def get_dels(wgdels_in):
    bed = []
    with open(wgdels_in, 'r') as dels:
        for l in dels:
            parts = l.rstrip().split('\t')
            bed.append(parts)

    return bed


def get_vars(wgvars_in):
    bed = []
    with open(wgvars_in, 'r') as vars:
        for l in vars:
            parts = l.rstrip().split('\t')
            bed.append([parts[4], parts[5], parts[7], parts[3], parts])

    return bed


def annotate_vars(dels, vars):
    b1 = get_dels(dels)
    b2 = get_vars(vars)

    annotated_vars = []

    for v in b2:
        c2, s2, e2, t, l = v
        for d in b1:
            c1, s1, e1 = d
            if t == 'DEL' and c1 == c2:
                if ( abs(int(s1) - int(s2)) < int(5e3) ) and ( abs(int(e1) - int(e2)) < int(5e3) ):
                    notes = l[26].split('; ')
                    notes.append('wg_del:True')
                    l[26] = '; '.join(notes)
            else:
                print("False negative: %s:%s%s in %s" % (c1, s1, e1, l[0]))
        annotated_vars.append('\t'.join(l))

    with open('data/all_WG_samples_merged_filt_annotated.txt', 'w') as f:
        f.write('\n'.join(map(str, annotated_vars)))


annotate_vars(wgdels_in, wgvars_in)