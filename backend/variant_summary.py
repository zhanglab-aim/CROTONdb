'''
Make table with genes with large difference
ZZJ
'''

import pandas as pd
import numpy as np


genenames = ['2/PDCD1', '2/CTLA4', '3/CCR5', 'X/ACE2']
#dfs = [pd.read_table("./data/data/genome/%s.tsv"%g).sort_values(by='abs_diff', ascending=False).head(n=3) for g in genenames]
dfs = [pd.read_table("./datavl/variant/tsv/%s.tsv"%g) for g in genenames]
for i in range(len(dfs)):
    dfs[i]['genename'] = genenames[i]
df = pd.concat(dfs, axis=0)

for x in ['ref_preds', 'alt_preds']:
    df[x] = df[x].apply(lambda x:
                               np.fromstring(
                               x.replace('\n','')
                                .replace('[','')
                                .replace(']','')
                                .replace('  ',' '), sep=' '))

tasks = ['Deletion Freq.', '1 bp Insertion', '1 bp Deleltion', '1 bp Frameshift Freq.',
        '2 bp Frameshift Freq.', 'Frameshift Freq.']


df_alt = df.loc[df['abs_diff']>0.25]
max_diff_tasks = []
ref_preds = []
alt_preds = []
for i in range(df_alt.shape[0]):
    diff = df_alt.iloc[i]['ref_preds'] - df_alt.iloc[i]['alt_preds']
    diff_idx = np.abs(diff).argmax()
    task = tasks[diff_idx]
    ref_preds.append(df_alt.iloc[i]['ref_preds'][diff_idx])
    alt_preds.append(df_alt.iloc[i]['alt_preds'][diff_idx])
    max_diff_tasks.append(task)

df_alt['max_diff_task'] = max_diff_tasks
df_alt['diff_ref_preds'] = ref_preds
df_alt['diff_alt_preds'] = alt_preds

df_alt_ = df_alt[['genename', 'chrom', 'start', 'end', 'strand', 'varid', 'pos', 'ref', 'alt', 'inx', 'max_diff_task', 'diff_ref_preds', 'diff_alt_preds', 'abs_diff', 'ref_preds', 'alt_preds']].sort_values(by='abs_diff', ascending=False)

df_alt_['Ref pred.'] = [','.join(['%.3f'%x for x in df_alt_.iloc[i]['ref_preds']])  for i in range(df_alt_.shape[0]) ]
df_alt_['Alt pred.'] = [','.join(['%.3f'%x for x in df_alt_.iloc[i]['alt_preds']])  for i in range(df_alt_.shape[0]) ]
df_alt_.drop(columns=['ref_preds', 'alt_preds'], inplace=True)

df_alt_.to_csv('gene_variant_summary.tsv', sep="\t", index=False)
