#!/usr/bin/env python
from __future__ import division
import os, re, sys
import json
from collections import defaultdict

test = False

# truth_set_in = 'data/visor_svs.bed'
truth_set_in = 'data/tumour_svs.bed'
vars_in = 'data/all_sim_samples_merged_reps_filt.txt'
if test:
    vars_in = 'data/R11_merged_filt_annotated.txt'

out_dir = 'data/'

#####################################
# Things that don't work:          ##
#   - no vars on Y...              ##
#####################################

fix_errors = True


def get_vars(simvars_in):
    bed = defaultdict(dict)
    svs = []

    with open(simvars_in, 'r') as vars:
        for l in vars:
            parts = l.rstrip().split('\t')
            orig_sv_type = parts[3]
            sv_type = parts[3]
            if parts[0] != 'sample':
                orig_sv_type = parts[-2].split(';')[0].split('=')[1]
            sample, c, s, e, t = parts[0], parts[4], parts[5], parts[7], convert_types(orig_sv_type)

            # if sample != 'visorR11': continue

            svs.append([c, s, e, t, parts])
            bed.setdefault(sample, []).append([c, s, e, t])

    write_bed(bed)

    return svs


def convert_types(sv):
    sv_types = {
        'BND': 'inversion',
        'INV': 'inversion',
        'DEL': 'deletion',
        'DUP': 'tandem duplication',
        'TANDUP': 'tandem duplication'
    }

    if sv in sv_types:
        return sv_types[sv]

    return sv.lower()


def get_truth(truth_file):

    d = defaultdict(lambda: defaultdict(dict))

    unmappable = {
        '2L_10625000': 1
    }

    count = 0
    with open(truth_file, 'r') as truth:
        for l in truth:
            parts = l.rstrip().split('\t')
            c, s, t = parts[0], parts[1], parts[3]
            k = '_'.join([c, s])
            if fix_errors:
                if c == 'Y': continue
                # if t == 'tandem duplication': continue
                if k in unmappable: continue

            d[c][s] = parts
            count += 1

    return d, count


def annotate_vars(vars, truth):
    vars = get_vars(vars)
    true_calls, true_call_count = get_truth(truth)

    annotated_vars = []
    true_positive_count = defaultdict(lambda: defaultdict(int))
    true_positive_types = defaultdict(lambda: defaultdict(int))
    false_positive_types = defaultdict(lambda: defaultdict(int))

    seen_events = defaultdict(lambda: defaultdict(int))
    keys = defaultdict(lambda: defaultdict(int))

    for v in vars:
        c, s, e, t, l = v
        sample, event, caller = l[:3]

        if sample == 'sample':
            annotated_vars.append('\t'.join(l))
            continue

        true_positive_count[sample]
        window = range(int(s) - 10, int(s) + 10, 1)

        if caller == 'Control-Freec':
            window = range(int(s) - 10000, int(s) + 10000, 1)

        for pos in iter(window):
            pos = str(pos)
            k = '_'.join([c, pos])
            true_positive, true_positive_count = is_true(true_calls, true_positive_count, c, pos, sample, t)
            if true_positive and not seen_events[sample][event]: break
        # keys['truth']

        if true_positive and not seen_events[sample][event]:
            keys[sample][k] = [event, 'tp', t]
            true_positive_types[sample][t] += 1

        elif not true_positive and not seen_events[sample][event]:
            keys[sample][k] = [event, 'fp', t]
            false_positive_types[sample][t] += 1

        seen_events[sample][event] += 1
        annotated_vars.append('\t'.join(l))

    keys = add_false_negatives(keys, true_calls)

    write_error_rates(keys)

    total_count = 0
    for s in true_positive_count:
        count = len(true_positive_count[s])
        total_count += count
        perc = (int(count)/true_call_count) * 100
        print("Sample: %s : %s / %s true positives found [%s]%%" % (s, count, true_call_count, int(perc)))

    no_samples = len(true_positive_count)
    print("Total: %s / %s [%s]%%" % (total_count, true_call_count*no_samples, round(total_count/(true_call_count*no_samples)*100)))

    with open('data/all_sim_samples_merged_filt_annotated.txt', 'w') as f:
        f.write('\n'.join(map(str, annotated_vars)))


def add_false_negatives(keys, true_calls):
    for s in keys:
        for c in true_calls:
            for pos in true_calls[c]:
                k = '_'.join([c, pos])
                # keys['truth'][k] = ['-', 'tp', true_calls[c][pos][3]]
                if k not in keys[s]:
                    # print("false negative: [%s] %s:%s-%s" % (true_calls[c][pos][3], true_calls[c][pos][0], true_calls[c][pos][1], true_calls[c][pos][2]))
                    keys[s][k] = ['-', 'fn', true_calls[c][pos][3]]

    return keys


def is_true(true_calls, true_positive_count, c, s, sample, sv_type):
    true_positive = False
    if c in true_calls and s in true_calls[c]:
        k = '_'.join([c, s])
        true_positive_count[sample][k] += 1
        true_positive = True
    return true_positive, true_positive_count


def write_bed(bed):
    for sample in bed.keys():
        if sample == 'sample': continue
        bed_out = os.path.join(out_dir, sample + '_svs.bed')
        with open(bed_out, 'w') as b:
            for c, s, e, t in bed[sample]:
                l = '\t'.join(map(str, [c, s, e, t]))
                b.write(l + '\n')
    return


def write_error_rates(vars):
    # print(json.dumps(vars, indent=4, sort_keys=True))

    df_out = os.path.join(out_dir, 'error_rates.txt')
    with open(df_out, 'w') as out:
        l = '\t'.join(['sample', 'event', 'chrom', 'start', 'type', 'error_group', 'rep', 'depth', 'purity'])
        out.write(l + '\n')

        sample_count = 0
        rep = 1
        purity_adj = 0

        for s in sorted(vars.keys(), key=lambda x: int(x.split('R')[1])):
            if sample_count % 5 == 0: purity_adj = 0

            sample_count += 1
            condition = 'shallow'
            if sample_count > 10:
                rep += 1
                sample_count = 1

            if sample_count > 5: condition = 'deep'

            purity = 100 - purity_adj
            purity_adj += 20

            print("[%s] %s: rep %s (%s). Purity: %s" % (sample_count, s, rep, condition, purity))
            for key in vars[s]:
                l = vars[s][key]
                e, et, sv = l
                c, pos = key.split('_')
                l = '\t'.join(map(str, [s, e, c, pos, sv, et, rep, condition, purity]))
                # print(l)
                out.write(l + '\n')
                # print("sample:%s, tp total: %s, tp types: %s, tp counts: %s, fp types: %s, fp counts: %s" % (s, len(tp_counts[s]), tp[s].keys(), tp[s].values(), fp[s].keys(), fp[s].values()))


annotate_vars(vars_in, truth_set_in)

