#!/usr/bin/env python
from __future__ import division
import os, re, sys
import json
from collections import defaultdict

wgdels_in = 'data/wholegut_dels_merged.bed'
truth_set_in = 'data/wholegut_dels_truth_set.bed'
wgvars_in = 'data/all_WG_samples_merged_filt.txt'
wgvars_in = 'data/all_WG_samples_merged_snvs_filt.txt'

def get_dels(dels_file):
    bed = []
    with open(dels_file, 'r') as dels:
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

def get_truth(truth_file):

    d = defaultdict(lambda: defaultdict(dict))

    with open(truth_file, 'r') as dels:
        for l in dels:
            parts = l.rstrip().split('\t')
            d[parts[0]][parts[1]] = parts

    return d

def annotate_vars(dels, vars, truth):
    dels = get_dels(dels)
    vars = get_vars(vars)
    true_calls = get_truth(truth)

    window = 2e5
    annotated_vars = []
    true_positive_count = defaultdict(lambda: defaultdict(int))

    for v in vars:
        true_positive = False
        wg_del = False
        c2, s2, e2, t, l = v
        if l[0] != 'sample': true_positive_count[l[0]]
        for d in dels:
            c1, s1, e1 = d
            if t == 'DEL' and c1 == c2:
                if ( abs(int(s1) - int(s2)) < window ) and ( abs(int(e1) - int(e2)) < window ):
                    notes = l[26].split('; ')
                    notes.append('wg_del:True')
                    l[26] = '; '.join(notes)
                    wg_del = True
                    if true_calls[c1] and true_calls[c1][s1]:
                        k = '_'.join([c1, s1])
                        true_positive_count[l[0]][k] += 1
                        true_positive = True

        if true_positive:
            print("True positive: [%s] %s:%s-%s", (l[0], c2, s2, e2))
        # if wg_del:
        #     print("Whole-gut deletion: [%s] %s:%s-%s", (l[0], c2, s2, e2))

        annotated_vars.append('\t'.join(l))

    print(json.dumps(true_positive_count, indent=2))

    total_count = 0
    for s in true_positive_count:
        count = len(true_positive_count[s])
        total_count += count
        perc = (int(count)/len(true_calls)) * 100
        print("Sample: %s : %s / %s true positives found [%s]%%" % (s, count, len(true_calls), perc))

    no_samples = len(true_positive_count)
    print("Total: %s / %s [%s]%%" % (total_count, len(true_calls)*no_samples, total_count/(len(true_calls)*no_samples)*100))

    with open('data/all_WG_samples_merged_filt_annotated.txt', 'w') as f:
        f.write('\n'.join(map(str, annotated_vars)))


annotate_vars(wgdels_in, wgvars_in, truth_set_in)