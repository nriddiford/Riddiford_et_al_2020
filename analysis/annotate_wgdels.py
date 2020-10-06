#!/usr/bin/env python
import os, re, sys
import ntpath
import pandas as pd
import csv

wgdels_in = 'data/wholegut_dels_merged.bed'
wgvars_in = 'data/all_WG_samples_merged_filt.txt'

def getDels(wgdels_in):
    bed = []
    with open(wgdels_in, 'r') as dels:
        for l in dels:
            parts = l.rstrip().split('\t')
            bed.append(parts)

    return bed


def getVars(wgvars_in):
    bed = []
    with open(wgvars_in, 'r') as vars:
        for l in vars:
            parts = l.rstrip().split('\t')
            bed.append([parts[4], parts[5], parts[7], parts[3], parts])

    return bed


def annotateDels(dels, vars):
    b1 = getDels(dels)
    b2 = getVars(vars)

    does_overlap = False

    for d in b1:
        c1, s1, e1 = d
        for v in b2:
            c2, s2, e2, t, l = v
            if t == 'DEL' and c1 == c2:
                if ( abs(int(s1) - int(s2)) < int(5e3) ) and ( abs(int(e1) - int(e2)) < int(5e3) ):
                    does_overlap = True


annotateDels(wgdels_in, wgvars_in)