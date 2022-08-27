"""get the genome-wide distribution statistics for a compiled tabix dir
"""

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import Counter, defaultdict
import pickle


TABIX_DIR = "./datavl/052121-variant-2/tabix/"
REF_PRED_COL = ['ref_1ins', 'ref_onemod3', 'ref_twomod3', 'ref_frameshift']
DIFF_PRED_COL = ['diff_1ins', 'diff_onemod3', 'diff_twomod3', 'diff_frameshift']

#CHRS = [str(i) for i in range(1, 23)] + ['X']
CHRS = ['22', '21']


def get_refpreds_varcnts(grouped_df):
    ref_preds = []
    var_cnts = []
    for _, grna_df in tqdm(grouped_df):
        ref_preds.append(['%.3f'%x for x in grna_df.iloc[0][REF_PRED_COL]])
        var_cnts.append(grna_df.shape[0])
    ref_preds = np.array(ref_preds)
    ref_preds = {x: Counter(ref_preds[:,i]) for i, x in enumerate(REF_PRED_COL)}
    var_cnts = Counter(var_cnts)
    return ref_preds, var_cnts


def get_diff_preds(chrom_df):
    diff_preds = [[]]*chrom_df.shape[0]
    for i, row in tqdm(chrom_df.iterrows(), total=chrom_df.shape[0]):
        #diff_preds.append(['%.3f'%abs(x) for x in row[DIFF_PRED_COL]])
        diff_preds[i] = ['%.3f'%abs(x) for x in row[DIFF_PRED_COL]]
    diff_preds = np.array(diff_preds)
    diff_preds = {x:Counter(diff_preds[:,i]) for i, x in enumerate(DIFF_PRED_COL)}
    return diff_preds


def merge_genome_wide_cdf(pred_dict):
    #merged_dict = {}
    running_sum = {}
    for i in np.arange(0, 1.001, 0.001):
        s = '%.3f' % i
        #merged_dict[s] = defaultdict(int)
        running_sum[s] = defaultdict(int)
        for chrom in CHRS:
            for item in pred_dict[chrom]:
                #merged_dict[s][item] += pred_dict[chrom][item][s]
                running_sum[s][item] += pred_dict[chrom][item][s]
        if i > 0:
            for item in running_sum[s]:
                running_sum[s][item] += running_sum['%.3f'% (i-0.001)][item]

    for i in np.arange(0, 1.001, 0.001):
        s = "%.3f" % i
        for item in running_sum[s]:
            running_sum[s][item] /= running_sum['1.000'][item]

    return running_sum


def get_per_bp_varcnt(var_cnts):
    tot = defaultdict(int)
    for chrom in CHRS:
        for varcnt in var_cnts[chrom]:
            tot[varcnt] += var_cnts[chrom][varcnt]
    normalizer = sum(tot.values())
    for varcnt in tot:
        tot[varcnt] /= normalizer
    avg = sum([k*v for k,v in tot.items()])
    per_bp = avg / 60
    return per_bp


def reload():
    gw_ref_preds = pickle.load(open('./data/052121-variant-2/gw_ref_preds.pkl', 'rb'))
    gw_diff_preds = pickle.load(open('./data/052121-variant-2/gw_diff_preds.pkl', 'rb'))
    gw_var_cnts = pickle.load(open('./data/052121-variant-2/gw_var_cnts.pkl', 'rb'))


def main():
    gw_ref_preds = {}
    gw_var_cnts = {}
    gw_diff_preds = {}
    for chrom in CHRS:
        print(chrom)
        dat = pd.read_table(os.path.join(TABIX_DIR, chrom, 'CROTON_varpred_%s.tsv'%chrom))
        diff_preds = get_diff_preds(dat)
        grouped_df = dat.groupby('pamid')
        ref_preds, var_cnts = get_refpreds_varcnts(grouped_df)
        # store
        gw_ref_preds[chrom] = ref_preds
        gw_var_cnts[chrom] = var_cnts
        gw_diff_preds[chrom] = diff_preds

    pickle.dump(gw_ref_preds, open('./data/052121-variant-2/gw_ref_preds.pkl', 'wb'))
    pickle.dump(gw_var_cnts, open('./data/052121-variant-2/gw_var_cnts.pkl', 'wb'))
    pickle.dump(gw_diff_preds, open('./data/052121-variant-2/gw_diff_preds.pkl', 'wb'))

    ref_cdf = merge_genome_wide_cdf(gw_ref_preds)
    diff_cdf = merge_genome_wide_cdf(gw_diff_preds)
    per_bp_varcnt_cdf = get_per_bp_varcnt(gw_var_cnts)

    data = {
        'ref_cdf': ref_cdf,
        'diff_cdf': diff_cdf,
        'perbp_varcnt_cdf': per_bp_varcnt_cdf
    }
    pickle.dump(data, open('./data/052121-variant-2/gw_data.pkl', 'wb'))

