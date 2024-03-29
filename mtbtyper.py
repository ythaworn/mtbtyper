#!/usr/bin/env python3
#
# mtbtyper
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Copyright (c) 2021
# Yuttapong Thawornwattana & Bharkbhoom Jaemsai
#
# Requirements:
# (1) python packages:
#     - numpy
#     - pandas
#     - scikit-allel
# (2) SNP schemes (csv files) in snpdb/ directory


import os
import re
import timeit
import argparse
import sys

import allel
import pandas as pd
import numpy as np


ref_lineage = ["L4", "L4.9", "L4.9(C)", "lineage4", "L4.2.1.1.1.1.1.1.i2"]
snp_freq_cutoff = 0.5  # cutoff for final prediction


def sort_by_freq(n_snp_list, n_snp_table):
    """Sort SNP counts by frequencies
    """

    # Convert series into dataframe
    n_snp_list_df = pd.DataFrame({'lineage': n_snp_list.index, 'n': n_snp_list.values})
    n_snp_table_df = pd.DataFrame({'lineage': n_snp_table.index, 'n_all': n_snp_table.values})
    
    n_snp = n_snp_list_df \
      .merge(n_snp_table_df, how='left', on='lineage') \
      .assign(freq=lambda x: x.n / x.n_all) \
      .sort_values(by=['n', 'freq'], ascending=False)
    
    return n_snp


def format_pred(pred):
    """Format list of genotype-specific SNPs as one string
    """

    if not pred.empty:
        pred_fmt = pred.lineage + ' (' + pred.n.astype(str) + '/' + pred.n_all.astype(str) + ')'
        pred_fmt = '"' + ', '.join(pred_fmt.tolist()) + '"'
    else:
        pred_fmt = 'unknown'

    return pred_fmt


def get_genotype_level(s):
    p_l2_modern = re.compile('2.2.M[1-6]') # higher than L2.2.Modern
    p_l2_2_1 = re.compile('2.2.1') # lower than L2.2.Modern

    out = s.count('.') + \
      (0.5 if re.search(p_l2_modern, s) else 0) - \
      (0.5 if re.search(p_l2_2_1, s) else 0) 
    return out


def predict_lineage_final(pred):
    """Make final lineage prediction
    """

    pred_final = 'unknown'

    if not pred.empty:
        pred_level = pred \
            .assign(level=[get_genotype_level(s) for s in pred.lineage]) \
            .sort_values(by='level', ascending=False)

        pred_level = pred_level[pred_level.freq > snp_freq_cutoff] \
            .reset_index(drop=True)

        if not pred_level.empty:
            pred_final = pred_level.lineage[0]

    return pred_final


def predict_lineage(snp_list, snp_table, n_snp_table, fotmat_output=False):
    """Find lineage-specific SNPs from a given SNP genotype scheme
    """

    pred = snp_list.merge(snp_table, on=['position', 'allele_change'])
    pred = pred[~pred.lineage.isin(ref_lineage)]  # exclude SNPs for ref lineage
    pred = pred[~pred.lineage.str.contains('*', regex=False)]  # for freschi2020
    out = pred.lineage.value_counts()

    # for ref lineage, use absense of SNP positions, ignoring allele changes
    pred_ref = snp_table.merge(snp_list, on='position', how='left', indicator=True)
    pred_ref = pred_ref[pred_ref.columns.drop(list(pred_ref.filter(regex='allele_change')))]
    pred_ref = pred_ref[pred_ref['_merge'] == 'left_only']
    ref_lineage_ind = pred_ref.lineage.isin(ref_lineage) | pred_ref['lineage'].str.contains('*', regex=False)
    pred_ref = pred_ref[ref_lineage_ind]

    if not pred_ref.empty:
        out_ref = pred_ref['lineage'].value_counts()
        out = pd.concat([out, out_ref])

    out = sort_by_freq(out, n_snp_table)
    
    if fotmat_output:
        out = format_pred(out)

    return out


def main(args):
    # check if input arguments are valid
    if not os.path.isdir(args.vcf_dir): 
        sys.exit('ERROR: invalid vcf dir: ' + args.vcf_dir)

    # output dir
    fout = os.path.join(args.out_dir, args.fout)
    os.makedirs(args.out_dir, exist_ok=True)

    # get list of vcf inputs
    f_vcf = [os.path.join(args.vcf_dir, f) for f in os.listdir(args.vcf_dir) if f.endswith(args.vcf_ending)]
    if not len(f_vcf) > 0: 
        sys.exit('no vcf files detected')
    
    # prepare SNP schemes
    snp_table = pd.read_csv(os.path.join(args.snpdb, 'main.csv'))
    n_snp_table = snp_table.lineage.value_counts()

    if args.all_schemes:
        all_schemes = [re.sub('\..*', '', f) for f in os.listdir(args.snpdb) if f.endswith('csv') and f != 'main.csv']
        all_schemes.sort()

        snp_tables = {}
        n_snp_tables = {}
        for i in range(len(all_schemes)):
            s = all_schemes[i]
            snp_tables[s] = pd.read_csv(os.path.join(args.snpdb, s + '.csv'))
            n_snp_tables[s] = snp_tables[s].lineage.value_counts()

    
    # write csv header line
    hdr = ','.join(['sample_id', 'genotype', 'genotype_specific_snp'])
    if args.all_schemes:
        hdr = hdr + ',' + ','.join(all_schemes)

    with open(fout, 'w') as file:
        file.write(hdr + '\n')

    # loop over vcf files
    for f in f_vcf:
        if not args.quiet: print(os.path.basename(f))

        # read SNPs from vcf
        f_tbi = f + '.tbi'
        callset = allel.read_vcf(f, numbers={'GT': 1}, tabix=f_tbi)
        h = allel.HaplotypeArray(callset['calldata/GT'])

        # exclude non-variant positions
        snp_ind = np.where(h.is_alt().ravel())[0]
        pos = callset['variants/POS'][snp_ind]
        ref = callset['variants/REF'][snp_ind]
        alt = callset['variants/ALT'][snp_ind, 0]

        snp_list = pd.DataFrame({'position': pos, 'allele_change': map('/'.join, zip(ref, alt))})

        # predict lineages
        snp_pred = predict_lineage(snp_list, snp_table, n_snp_table)
        snp_pred_fmt = format_pred(snp_pred)

        # make final prediction
        final_genotype = predict_lineage_final(snp_pred)

        if args.all_schemes:
            snp_pred_all = [None] * len(all_schemes)
            for i in range(len(all_schemes)):
                s = all_schemes[i]
                pred = predict_lineage(snp_list, snp_tables[s], n_snp_tables[s], fotmat_output=True)
                snp_pred_all[i] = pred

        # output
        sample_id = os.path.basename(f)
        sample_id = re.sub('\..*', '', sample_id)

        with open(fout, 'a') as file:
            if args.all_schemes:
                file.write(','.join([sample_id, final_genotype, snp_pred_fmt] + snp_pred_all) + '\n')
            else:
                file.write(','.join([sample_id, final_genotype, snp_pred_fmt]) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Predict Mtb lineage.')

    # Required positional argument
    parser.add_argument('vcf_dir', metavar='vcf_dir', type=str,
                        help='directory to vcf files')

    # Optional arguments
    parser.add_argument('-o', '--out', type=str, 
        dest='out_dir', default=os.getcwd(),
        help='output directory (default: current working directory)')

    parser.add_argument('-f', '--fout', type=str, 
        dest='fout', default='lineage.csv',
        help='output file name (default: lineage.csv)')

    parser.add_argument('-e', '--vcf_end', type=str, 
        dest='vcf_ending', default='vcf.gz',
        help='ending pattern of vcf file (default: vcf.gz)')

    parser.add_argument('--all_schemes', action='store_true',
        help='Add prediction from all available SNP schemes (default: false)')

    parser.add_argument('--snpdb', type=str, 
        dest='snpdb', default='snpdb',
        help='Path to genotyping SNP schemes (default: snpdb/)')

    parser.add_argument('--quiet', action='store_true',
        help='Suppress screen output (default: false)')

    args = parser.parse_args()

    start = timeit.default_timer()
    main(args)
    elapsed = timeit.default_timer() - start
    if not args.quiet: print("elapsed time: %.2f s" % elapsed)
