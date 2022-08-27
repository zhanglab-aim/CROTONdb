'''
Compare database versions
041321
'''

import os
import pandas as pd


bed_dir = './datavl/variant/tsv'
#bed_dir_og = './datavl/variant-before/tsv'
fp_lst = sorted([os.path.join(path, name) for path, subdirs, files in os.walk(bed_dir) for name in files])
different = []
for fp in fp_lst:
    df = pd.read_csv(fp, sep='\t')
    fp_og = fp.replace('/variant/', '/variant-before/')
    df_og = pd.read_csv(fp_og, sep='\t')

    if not df.equals(df_og): different.append(fp)
       
print(different)

different = ['./datavl/variant/tsv/5/ADGRV1.tsv', './datavl/variant/tsv/5/AMACR.tsv', './datavl/variant/tsv/5/ANXA2R.tsv', './datavl/variant/tsv/5/ARL10.tsv', './datavl/variant/tsv/6/ECT2L.tsv', './datavl/variant/tsv/6/NDUFAF4.tsv', './datavl/variant/tsv/6/SOD2.tsv', './datavl/variant/tsv/7/ABCF2-H2BE1.tsv', './datavl/variant/tsv/7/COBL.tsv', './datavl/variant/tsv/7/GPR85.tsv', './datavl/variant/tsv/7/NCAPG2.tsv', './datavl/variant/tsv/7/SAMD9L.tsv', './datavl/variant/tsv/7/TRBJ2-6.tsv', './datavl/variant/tsv/7/ZNF746.tsv', './datavl/variant/tsv/8/DEFA5.tsv', './datavl/variant/tsv/8/KIFC2.tsv', './datavl/variant/tsv/8/PSKH2.tsv', './datavl/variant/tsv/8/UBXN2B.tsv', './datavl/variant/tsv/9/C9orf153.tsv', './datavl/variant/tsv/9/FOXD4L4.tsv', './datavl/variant/tsv/9/NDUFA8.tsv', './datavl/variant/tsv/9/RORB.tsv', './datavl/variant/tsv/9/UNC13B.tsv', './datavl/variant/tsv/X/CHM.tsv', './datavl/variant/tsv/X/GAGE2A.tsv', './datavl/variant/tsv/X/MOSPD2.tsv', './datavl/variant/tsv/X/RPL39.tsv', './datavl/variant/tsv/X/UBE2A.tsv']

