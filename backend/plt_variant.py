'''
Functions to plot figs for ISMB/ECCB submission
012721
'''

import sys
#sys.path.insert(0,'/Users/victoriali/Documents/GitHub/CROTON')
#sys.path.insert(0,'/mnt/ceph/users/zzhang/croton-2')
sys.path.insert(0,'/mnt/ceph/users/vli/croton')
from src.variant.get_CDSpams import get_CDS_for_gene
from src.variant.make_tsvs import get_varscore
from src.variant.make_tsvs import filter_disruptive_pams
import matplotlib.pyplot as plt
import seaborn as sns
import gffutils
import numpy as np
import pandas as pd

import bisect
import intervaltree

# ALSO SEE images/ISMB_ECCB/revision.py for the PDCD1 plot in the manuscript
#croton1 GACTGCCGCTTCCGTGTCACACAACTGCCCAACGGGCGTGACTTCCACATGAGCGTGGTC
#croton2                       AACTGCCCAACGGGCGTGACTTCCACATGAGCGTGGTCAGGGCCCGGCGCAATGACAGCG

db_path = '/mnt/ceph/users/zzhang/decan/data/genomes/gencode/gencode.v35.annotation.gff3.gz.decan_db'
db = gffutils.FeatureDB(db_path)
gn_dict = {g.attributes['gene_name'][0]: g.id for g in db.features_of_type('gene')}

def insert_exonic_region_number(df, cds_array, strand):
    exonic_region_num = []
    for i in range(df.shape[0]):
        ctr = (df.iloc[i]['start'] + df.iloc[i]['end'] ) //2
        ctr = ctr + 5 if strand =='+' else ctr - 5
        exonic_region_num.append(bisect.bisect(cds_array, ctr))
    df['exonic_region'] = exonic_region_num
    return df

def dotplt(genename):
    global db, gn_dict
    bed_df = filter_disruptive_pams(genename)
    cds_array = get_CDS_for_gene(genename)
    gene = db[gn_dict[genename]]
    bed_df = insert_exonic_region_number(bed_df, cds_array, strand=gene.strand)
    bed_df = bed_df.sort_values(by=['start'], ascending=gene.strand=='+')
    if gene.strand=='-':
        bed_df['exonic_region'] = len(cds_array) - bed_df['exonic_region']
    bed_df['Exonic Region'] = ['E%i'%i for i in bed_df['exonic_region']]

    # this will create boxes for PAMs on both +/- strands
    #start_lst = list(bed_df['start'].unique())
    #for i in range(len(start_lst)):
    bed_df.index = pd.MultiIndex.from_tuples(list(zip(bed_df['start'], bed_df['strand'])))
    bed_df['pam_num'] = np.nan
    for i, idx in enumerate(bed_df.index.unique()):
        #bed_df.at[bed_df.start == start_lst[i], 'pam_num'] = i + 1
        bed_df.at[idx, 'pam_num'] = i + 1
    bed_df['pam_num'] = bed_df['pam_num'].astype(int)
    ref_df = bed_df[['pam_num', 'ref_frameshift']]
    ref_df = ref_df.rename(columns={'ref_frameshift':'Frameshift Frequency'})
    ref_df = ref_df.drop_duplicates()
    ref_df['type'] = 'reference'
    alt_df = bed_df[['pam_num', 'alt_frameshift', 'Exonic Region']]
    alt_df = alt_df.rename(columns={'alt_frameshift':'Frameshift Frequency'})
    alt_df['type'] = 'alternate'
    plt_df = pd.concat([ref_df, alt_df])
 
    plt.figure(figsize=(24,6))
    sns.set(style="whitegrid")
    sns.stripplot(x='pam_num', y='Frameshift Frequency', hue='Exonic Region', data=alt_df, marker='o', linewidth=1) #, jitter=False, dodge=True, , color='black'
    ax = sns.boxplot(x='pam_num', y='Frameshift Frequency', data=ref_df, showfliers=False, linewidth=3, color='black') #, palette=sns.color_palette(colors)
    ax.set_xlabel('PAM index')
    for spine in ['right', 'top', 'left', 'bottom']: ax.spines[spine].set_visible(False)
    plt.savefig('./images/ISMB_ECCB/gene_variants/%s-genewide.pdf' %(genename), bbox_inches='tight')#linewidth=1
    return plt_df

MEDIUM_SIZE = 12

def dotplt_delfreq(genename):
    global db, gn_dict
    bed_df = filter_disruptive_pams(genename)
    cds_array = get_CDS_for_gene(genename)
    gene = db[gn_dict[genename]]
    bed_df = insert_exonic_region_number(bed_df, cds_array, strand=gene.strand)
    bed_df = bed_df.sort_values(by=['start'], ascending=gene.strand=='+')
    if gene.strand=='-':
        bed_df['exonic_region'] = len(cds_array) - bed_df['exonic_region']
    bed_df['Exonic Region'] = ['E%i'%i for i in bed_df['exonic_region']]
    
    bed_df.index = pd.MultiIndex.from_tuples(list(zip(bed_df['start'], bed_df['strand'])))
    bed_df['pam_num'] = np.nan
    for i, idx in enumerate(bed_df.index.unique()):
        bed_df.at[idx, 'pam_num'] = i + 1
    bed_df['pam_num'] = bed_df['pam_num'].astype(int)
    print(bed_df[bed_df['ref_seq'] == 'AACTGCCCAACGGGCGTGACTTCCACATGAGCGTGGTCAGGGCCCGGCGCAATGACAGCG'])
    #GACTGCCGCTTCCGTGTCACACAACTGCCCAACGGGCGTGACTTCCACATGAGCGTGGTC
    ref_df = bed_df[['pam_num', 'ref_del_frq']]
    ref_df = ref_df.rename(columns={'ref_del_frq':'Deletion Frequency'})
    ref_df = ref_df.drop_duplicates()
    ref_df['type'] = 'reference'
    alt_df = bed_df[['pam_num', 'alt_del_frq', 'Exonic Region']]
    alt_df = alt_df.rename(columns={'alt_del_frq':'Deletion Frequency'})
    alt_df['type'] = 'alternate'
    plt_df = pd.concat([ref_df, alt_df])
    
    plt.figure(figsize=(24,6))
    sns.set(style="whitegrid")
    sns.boxplot(x='pam_num', y='Deletion Frequency', data=ref_df, showfliers=False, linewidth=3, color='black')
    ax = sns.stripplot(x='pam_num', y='Deletion Frequency', hue='Exonic Region', data=alt_df, marker='o', linewidth=1) 
   
    ax.set_xlabel('PAM index', fontsize=MEDIUM_SIZE)
    for spine in ['right', 'top', 'left', 'bottom']: ax.spines[spine].set_visible(False)
    plt.ylabel('Deletion Frequency', fontsize=MEDIUM_SIZE)
    plt.tick_params(labelsize=MEDIUM_SIZE-2)
    plt.legend(title='Exonic Region', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=MEDIUM_SIZE)
    plt.savefig('./images/ISMB_ECCB/horiz-%s-delfreq-genewide.png' %(genename), bbox_inches='tight', dpi=350)#linewidth=1
    return plt_df

LARGE_SIZE = 24
def boxplt_delfreq():
    task = 'del_frq'
    model_fp = "./outputs/forecast_freqs-4/train.bak6.wsf6/bestmodel.h5"
    df = get_varscore(bed_fp="./data/data/genome/PDCD1.bed", model_fp = model_fp)
    plt.figure(figsize=(15,8))
    ax = sns.boxplot(data=df, x='inx', y=df['diff_%s'%task].abs(), showfliers=False)
    plt.tick_params(labelsize=LARGE_SIZE-2)
    plt.xlabel('Position', fontsize=LARGE_SIZE)
    plt.xticks(ticks = np.append(np.arange(0, 60, step=5), 60), labels = np.append(np.arange(0, 60, step=5), 60))
    plt.ylabel('Deletion Frequency Difference', fontsize=LARGE_SIZE)
    plt.savefig('./images/ISMB_ECCB/gene_variants/%s_absdiff_box_%s.png'%(task, 'PDCD1'), dpi=350, bbox_inches='tight')

MEDIUM_SIZE = 20
def dotplt_1bpins(genename):
    global db, gn_dict
    bed_df = filter_disruptive_pams(genename)
    cds_array = get_CDS_for_gene(genename)
    gene = db[gn_dict[genename]]
    bed_df = insert_exonic_region_number(bed_df, cds_array, strand=gene.strand)
    bed_df = bed_df.sort_values(by=['start'], ascending=gene.strand=='+')
    if gene.strand=='-':
        bed_df['exonic_region'] = len(cds_array) - bed_df['exonic_region']
    bed_df['Exonic Region'] = ['E%i'%i for i in bed_df['exonic_region']]
    
    bed_df.index = pd.MultiIndex.from_tuples(list(zip(bed_df['start'], bed_df['strand'])))
    bed_df['pam_num'] = np.nan
    for i, idx in enumerate(bed_df.index.unique()):
        bed_df.at[idx, 'pam_num'] = i + 1
    bed_df['pam_num'] = bed_df['pam_num'].astype(int)

    ref_df = bed_df[['pam_num', 'ref_1ins']]
    ref_df = ref_df.rename(columns={'ref_1ins':'1 bp Insertion Probability'})
    ref_df = ref_df.drop_duplicates()
    ref_df['type'] = 'reference'
    alt_df = bed_df[['pam_num', 'alt_1ins', 'Exonic Region']]
    alt_df = alt_df.rename(columns={'alt_1ins':'1 bp Insertion Probability'})
    alt_df['type'] = 'alternate'
    plt_df = pd.concat([ref_df, alt_df])
    
    plt.figure(figsize=(24,6))
    sns.set(style="whitegrid")
    sns.boxplot(x='pam_num', y='1 bp Insertion Probability', data=ref_df, showfliers=False, linewidth=3, color='black')
    ax = sns.stripplot(x='pam_num', y='1 bp Insertion Probability', hue='Exonic Region', data=alt_df, marker='o', linewidth=0.5, size=6) 
    
    ax.set_xlabel('PAM index', fontsize=MEDIUM_SIZE)
    for spine in ['right', 'top', 'left', 'bottom']: ax.spines[spine].set_visible(False)
    plt.ylabel('1 bp Insertion Probability', fontsize=MEDIUM_SIZE)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.tick_params(labelsize=MEDIUM_SIZE-2)
    legend = plt.legend(title='Exonic\nRegion', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=MEDIUM_SIZE, handletextpad=0.1)
    plt.setp(legend.get_title(),fontsize=MEDIUM_SIZE)
    plt.savefig('./images/ISMB_ECCB/1bpins-genewide.png', bbox_inches='tight', dpi=350)#linewidth=1
    return plt_df

def plt_varscore(df, plttype, genename, task='frameshift'):
    plt.close()
    plt.clf()
    if plttype == 'scatter': # ref vs alt
        plt.figure(figsize=(8,8))
        plt.grid(True)
        ax = sns.scatterplot(data=df, x='ref_%s'%task, y='alt_%s'%task, markers=['.'])
        _min, _max = min(df['ref_%s'%task].min(), df['alt_%s'%task].min()), \
                    max(df['ref_%s'%task].max(), df['alt_%s'%task].max()),
        ax.plot([_min, _max], [_min, _max], linestyle='--', color='grey')
        plt.savefig('./images/ISMB_ECCB/gene_variants/scatter_%s_%s.png'%(task, genename))
    # score difference over position on sequence
    elif plttype == 'bar':
        avg_absdiff = []
        inxs = list(set(df['inx'].tolist()))
        for inx in inxs:
            df_inx = df[df['inx'] == inx]
            total_absdiff = df_inx['diff_%s'%task].abs().to_numpy().sum()
            avg_absdiff.append(total_absdiff / len(df_inx))
        inx_score_dict = {'inx':inxs, 'avg_absdiff':avg_absdiff}
        inx_score_df = pd.DataFrame(inx_score_dict)
        plt.figure(figsize=(25,10))
        sns.barplot(data=inx_score_df, x='inx', y='avg_absdiff')
        plt.savefig('./images/ISMB_ECCB/gene_variants/%s_absdiff_%s.png'%(task, genename))
    elif plttype == 'box':
        plt.figure(figsize=(14,8))
        ax = sns.boxplot(data=df, x='inx', y=df['diff_%s'%task].abs(), showfliers=False)
        plt.savefig('./images/ISMB_ECCB/gene_variants/%s_absdiff_box_%s.png'%(task, genename))
    elif plttype == 'violin':
        plt.figure(figsize=(14,8))
        ax = sns.violinplot(data=df, x='inx', y=df['diff_%s'%task].abs(), showfliers=False)
        plt.savefig('./images/ISMB_ECCB/gene_variants/%s_absdiff_violin_%s.png'%(task, genename))

def main():
    model_fp = "./outputs/forecast_freqs-4/train.bak6.wsf6/bestmodel.h5"
    genenames = ['PDCD1', 'CTLA4', 'CCR5', 'ACE2']
    plttypes = ['scatter', 'box']
    tasks = ['1ins', '1del', 'frameshift']
    for genename in genenames:
        print(genename)
        df = get_varscore(bed_fp="./data/data/genome/%s.bed"%genename,
                model_fp = model_fp)
        for plttype in plttypes:
            for task in tasks:
                plt_varscore(df, plttype, genename, task=task)
        df.to_csv('./data/data/genome/%s.tsv'%genename, index=False, sep="\t")
        df['genename'] = genename

    dfs = [filter_disruptive_pams(bed_df=pd.read_table('./data/data/genome/%s.tsv'%genename)) for genename in genenames]
    for i, df in enumerate(dfs):
        df['Gene'] = genenames[i]
    tot_df = pd.concat(dfs, axis=0)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    sns.scatterplot(data=tot_df, x='ref_frameshift', y='alt_frameshift', hue='Gene',
            alpha=0.75,
            hue_order=genenames,
            ax=ax)
    ax.grid(True)
    ax.plot([0.3,1], [0.3,1], linestyle='--', color='black')
