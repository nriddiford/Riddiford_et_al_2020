#!/usr/bin/env python
from __future__ import division
import os, re, sys
import json
from collections import defaultdict


truth_set_in = 'data/visor_svs.bed'
vars_in = 'data/all_sim_samples_merged_filt.txt'


def get_vars(simvars_in):
    bed = []
    with open(simvars_in, 'r') as vars:
        for l in vars:
            parts = l.rstrip().split('\t')
            bed.append([parts[4], parts[5], parts[7], parts[3], parts])

    return bed


def get_truth(truth_file):

    d = defaultdict(lambda: defaultdict(dict))

    count = 0
    with open(truth_file, 'r') as truth:
        for l in truth:
            count += 1
            parts = l.rstrip().split('\t')
            d[parts[0]][parts[1]] = parts

    return d, count


def annotate_vars(vars, truth):
    vars = get_vars(vars)
    true_calls, true_call_count = get_truth(truth)

    annotated_vars = []
    true_positive_count = defaultdict(lambda: defaultdict(int))
    seen_events = defaultdict(lambda: defaultdict(int))

    for v in vars:
        c, s, e, t, l = v
        sample, event, caller = l[:3]

        if sample == 'sample':
            annotated_vars.append('\t'.join(l))
            continue

        true_positive_count[sample]
        window = range(int(s) - 10000, int(s) + 10000, 1)

        for pos in window:
            pos = str(pos)
            true_positive = False
            if true_calls[c] and true_calls[c][pos]:
                k = '_'.join([c, pos])
                true_positive_count[sample][k] += 1
                true_positive = True

            # if true_positive:
            #     print("True positive: [%s] %s:%s-%s", (sample, c, pos, e))

        annotated_vars.append('\t'.join(l))

    # print(json.dumps(true_positive_count, indent=2))

    total_count = 0
    for s in true_positive_count:
        count = len(true_positive_count[s])
        total_count += count
        perc = (int(count)/true_call_count) * 100
        print("Sample: %s : %s / %s true positives found [%s]%%" % (s, count, true_call_count, perc))

    no_samples = len(true_positive_count)
    print("Total: %s / %s [%s]%%" % (total_count, true_call_count*no_samples, total_count/(true_call_count*no_samples)*100))

    with open('data/all_sim_samples_merged_filt_annotated.txt', 'w') as f:
        f.write('\n'.join(map(str, annotated_vars)))


annotate_vars(vars_in, truth_set_in)

