"""plot functions for vis in crotondb
FZZ
9/30/2021
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_AF(ele, return_df=False):
    tmp = [x.replace("'", '').replace('{','').replace('}','') for x in ele.split(",") ]
    info = {x.split(':')[0].strip(): float(x.split(':')[1]) for x in tmp}
    if return_df is True:
        info = pd.DataFrame.from_dict(info, orient='index', columns=['MAF'])
    return info


def get_candidates(df):
    cands = []
    for i, row in df.iterrows():
        if row['ref_frameshift'] > 0.9 or row['alt_frameshift'] > 0.9:
            cands.append(i)
        #af_info = parse_AF(row['AF_info'])
        #if not 'AF_eas' in af_info:
        #    continue
        #if af_info['AF_eas'] > 0.05:
        #    cands.append(i)
    return cands



def plot(pamid, df):
    row = df.query('pamid=="%s"'%pamid)
    score_cols = [f'{t}_{s}' for t in ['ref','alt'] for s in 
            #['1ins', '1del', 'onemod3', 'twomod3', 'frameshift']
            ['1ins', '1del', 'frameshift']
            ]
    score_df = row[score_cols].melt()
    score_df['allele'] = [x.split('_')[0] for x in score_df['variable']]
    score_df['allele'].replace({
        'ref': f'{row.iloc[0]["ref"]}(ref)',
        'alt': f'{row.iloc[0]["alt"]}(alt)',
        }, inplace=True
    )
    score_df['task'] = [x.split('_')[1] for x in score_df['variable']]

    af_df = parse_AF(row.iloc[0]['AF_info'], return_df=True)
    af_df = af_df.loc[
        ['AF_afr', 'AF_amr', 'AF_asj', 'AF_eas', 'AF_sas', 'AF_fin', 'AF_nfe', 'AF_oth']
    ].reset_index()

    # better naming
    score_df['task'].replace({
        '1del': '1bp Del',
        '1ins': '1bp Ins',
        'onemod3': '1bp Frameshift',
        'twomod3': '2bp Frameshift',
        'frameshift': 'Total Frameshift'
        }, inplace=True)
    af_df['index'].replace({
        'AF_afr': 'African/African American',
        'AF_amr': 'American Admixed/Latino',
        'AF_asj': 'Ashkenazi Jewish',
        'AF_eas': 'East Asian',
        'AF_fin': 'Finnish',
        'AF_nfe': 'Non-Finnish European',
        'AF_oth': 'Other',
        'AF_sas': 'South Asian'
        }, inplace=True)

    fig, axs = plt.subplots(1, 2, figsize=(10,6))
    sns.barplot(x='task', y='value', hue='allele', data=score_df, ax=axs[0])
    axs[0].set_xlabel('Editing Outcomes', fontsize='x-large', fontweight='light')
    axs[0].set_ylabel('Predicted Probability', fontsize='x-large', fontweight='light')
    axs[0].set_xticklabels(axs[0].get_xticklabels(),rotation = 30, horizontalalignment='right', fontsize='x-large', fontweight='light')
    axs[0].set_ylim(0, 1)
    sns.barplot(x='index', y='MAF', color="darkgray", data=af_df, ax=axs[1])
    axs[1].set_xlabel('Population', fontsize='x-large', fontweight='light')
    axs[1].set_ylabel('Minor Allele Frequency', fontsize='x-large', fontweight='light')
    axs[1].set_xticklabels(axs[1].get_xticklabels(),rotation = 30, horizontalalignment='right', fontsize='x-large', fontweight='light')
    fig.suptitle('{gene}\n{rs}, PAM {id}'.format(gene=pamid.split('|')[0], id=pamid.split('|')[1],
        rs=row.iloc[0]['id']))
    fig.tight_layout(rect=[0, 0.03, 1, 0.9])
    fig.savefig(f'{pamid}.png')


df = pd.read_table('./datavl/variant/filt-fs-AF.fs/variable-AF/merged.tsv')
pamids = [
    "TNFRSF25|160",
    "MUC16|2191",
    "FGFR3|147"]
for pamid in pamids:
    plot(pamid, df)


