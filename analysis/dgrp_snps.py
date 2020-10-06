#!/usr/bin/env python
import os, re, sys
import ntpath
import pandas as pd
import csv

snps_in='dm_DGRP_SNP.sga.gz'

def getSNPS(snps_in):
    df = pd.read_csv(snps_in, delimiter="\t", compression='gzip', index_col=False, na_filter=False, names=['id', 'type', 'end', 'strand', 'something', 'locus'])
    df[['chrom', 'start', 'detail1', 'info2']] = df['locus'].str.split('_', expand=True)

    chroms = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
    df['chrom'].isin(chroms)


    df = df[df.type == "snp"]

    df['end'] = df['start'].astype(int) + 1


    selection = df[['chrom', 'start', 'end']]


    # print(selection.describe(include='all'))

    selection.to_csv('dgrp_snps.bed', sep="\t", index=False)


    # with open('dgrp_snps.bed', 'w') as out_file:
    #     out_file.write('\t'.join(['chrom', 'start', 'end']) + '\n')
    #     for index, row in df.iterrows():
    #         chrom = row['chrom']
    #         start = int(row['start'])
    #         end = int(row['start']) + 1
    #
    #         line = '\t'.join(map(str, [chrom, start, end]))
    #         out_file.write(line + '\n')

getSNPS(snps_in)