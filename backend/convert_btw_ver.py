'''
Convert new tabix files to be similar to old version
VRL 080621
'''

import pandas as pd

# Load beginning of old_df
old_chr21_fp = 'datavl/052121-variant-2/tabix/21/CROTON_varpred_21.tsv' # Old version of df server is based on
old_chr21_df = pd.read_csv(old_chr21_fp, sep='\t', nrows=10)

# Alter new_df to include all rows of old_df
new_chr21_fp = 'datavl/variant/tabix/21/CROTON_varpred_21.tsv'
new_chr21_df = pd.read_csv(new_chr21_fp, sep='\t')

new_chr21_df[['genename','num']] = new_chr21_df.pamid.str.split('|', expand=True)
new_chr21_df['inx'] = new_chr21_df['pos'] - new_chr21_df['start']

tasks_ordering = ['del_frq', '1ins', '1del', 'onemod3', 'twomod3', 'frameshift']
for task in tasks_ordering:
    new_chr21_df['diff_%s'%task] = new_chr21_df['ref_%s'%task] - new_chr21_df['alt_%s'%task]

new_chr21_df['ref_preds'] = new_chr21_df[['ref_' + task for task in tasks_ordering]].values.tolist()
new_chr21_df['alt_preds'] = new_chr21_df[['alt_' + task for task in tasks_ordering]].values.tolist()

# Reorder columns in new_df by old_df order
old_order = old_chr21_df.columns.tolist()
old_order[0] = '#chrom'
old_order.append('AF_info') # Include new allele freq information column
new_chr21_df = new_chr21_df[old_order]

# Save file 
# NOTE: try no hashtag to see if that works
print(new_chr21_df)
new_chr21_df.to_csv('datavl/variant/prev-ver/21/CROTON_varpred_21.tsv', sep='\t', index=False)

'''
bgzip -c datavl/variant/prev-ver/21/CROTON_varpred_21.tsv > datavl/variant/prev-ver/21/CROTON_varpred_21.gz
tabix -p vcf datavl/variant/prev-ver/21/CROTON_varpred_21.gz
'''

#new_chr21_df = new_chr21_df.rename(columns={"#chrom": "chrom"}) 
# NOTE: HASHTAG IS NECESSARY in chrom column for tabix
    # BUT not included in prev version? that's confusing
'''
print(new_chr21_df.columns)
Index(['#chrom', 'pos', 'start', 'end', 'strand', 'pamid', 'id', 'ref', 'alt',
       'AF_info', 'ref_seq', 'alt_seq', 'abs_diff', 'ref_del_frq',
       'alt_del_frq', 'ref_1ins', 'alt_1ins', 'ref_1del', 'alt_1del',
       'ref_onemod3', 'alt_onemod3', 'ref_twomod3', 'alt_twomod3',
       'ref_frameshift', 'alt_frameshift'],
      dtype='object')

print(old_chr21_df.columns)
Index(['chrom', 'pos', 'id', 'ref', 'alt', 'pamid', 'genename', 'num', 'start',
       'end', 'strand', 'inx', 'ref_seq', 'alt_seq', 'ref_preds', 'alt_preds',
       'abs_diff', 'ref_del_frq', 'alt_del_frq', 'diff_del_frq', 'ref_1ins',
       'alt_1ins', 'diff_1ins', 'ref_1del', 'alt_1del', 'diff_1del',
       'ref_onemod3', 'alt_onemod3', 'diff_onemod3', 'ref_twomod3',
       'alt_twomod3', 'diff_twomod3', 'ref_frameshift', 'alt_frameshift',
       'diff_frameshift'],
      dtype='object')
'''
