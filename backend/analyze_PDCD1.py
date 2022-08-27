'''
Analyze PDCD1 for wet lab targeting

'''

import sys
#sys.path.insert(0,'/Users/victoriali/Documents/GitHub/CROTON')
#sys.path.insert(0,'/mnt/ceph/users/zzhang/croton-2')
sys.path.insert(0,'/mnt/ceph/users/vli/croton')

import pandas as pd
from src.variant.get_CDSpams import get_CDS_for_gene
from src.variant.plt_variant import insert_exonic_region_number
from src.variant.make_tsvs import get_varscore #filter_disruptive_pams
import gffutils
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pyfaidx import Fasta


db_path = '/mnt/ceph/users/zzhang/decan/data/genomes/gencode/gencode.v35.annotation.gff3.gz.decan_db'
MEDIUM_SIZE = 10
genome_path = './data/data/genome/GRCh38.primary_assembly.genome.fa'

def get_beddf():
    db = gffutils.FeatureDB(db_path)
    gn_dict = {g.attributes['gene_name'][0]: g.id for g in db.features_of_type('gene')}
    bed_df = pd.read_csv('datavl/variant/tsv/2/PDCD1.tsv', sep='\t')
    posdrop_series = (bed_df['strand']=='+') & ((bed_df['inx'] == 35)|(bed_df['inx'] == 36)) # interferes with PAM on + strand seqs
    bed_df = bed_df[~posdrop_series]
    negdrop_series = (bed_df['strand']=='-') & ((bed_df['inx'] == 25)|(bed_df['inx'] == 26)) # interferes with PAM on - strand seqs
    bed_df = bed_df[~negdrop_series]
    cds_array = get_CDS_for_gene('PDCD1')
    gene = db[gn_dict['PDCD1']]
    bed_df = insert_exonic_region_number(bed_df, cds_array, strand=gene.strand)
    bed_df = bed_df.sort_values(by=['start'], ascending=gene.strand=='+')
    if gene.strand=='-':
        bed_df['exonic_region'] = len(cds_array) - bed_df['exonic_region']

    #bed_df.index = pd.MultiIndex.from_tuples(list(zip(bed_df['start'], bed_df['strand'])))
    #bed_df['pam_num'] = np.nan
    #for i, idx in enumerate(bed_df.index.unique()):
    #    bed_df.at[idx, 'pam_num'] = i + 1
    #bed_df['pam_num'] = bed_df['pam_num'].astype(int)
    #bed_df = bed_df[['strand', 'pam_num', 'pos', 'varid', 'ref_seq', 'alt_seq', 
    #       'ref_del_frq', 'alt_del_frq', 'diff_del_frq', 'ref_1ins', 'alt_1ins', 'diff_1ins']]
    #bed_df.to_csv('./data/PDCD1.csv', index=False)

    return bed_df

def dotplt():
    bed_df = get_beddf()
    bed_df = bed_df[bed_df['exonic_region'] == 2]
    bed_df['Exonic Region'] = ['E%i'%i for i in bed_df['exonic_region']]
    ref_df = bed_df[['num', 'ref_del_frq']]
    ref_df = ref_df.rename(columns={'ref_del_frq':'Deletion Frequency'})
    ref_df = ref_df.drop_duplicates()
    ref_df['type'] = 'reference'
    alt_df = bed_df[['num', 'alt_del_frq', 'Exonic Region']]
    alt_df = alt_df.rename(columns={'alt_del_frq':'Deletion Frequency'})
    alt_df['type'] = 'alternate'
    plt_df = pd.concat([ref_df, alt_df])

    plt.figure(figsize=(24,6))
    sns.set(style="whitegrid")
    sns.boxplot(x='num', y='Deletion Frequency', data=ref_df, showfliers=False, linewidth=3, color='black')
    ax = sns.stripplot(x='num', y='Deletion Frequency', hue='Exonic Region', data=alt_df, marker='o', linewidth=0.5, size=6) 

    ax.set_xlabel('PAM index', fontsize=MEDIUM_SIZE)
    for spine in ['right', 'top', 'left', 'bottom']: ax.spines[spine].set_visible(False)
    plt.ylabel('Deletion Frequency', fontsize=MEDIUM_SIZE)
    plt.tick_params(labelsize=MEDIUM_SIZE-2)
    legend = plt.legend(title='Exonic\nRegion', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=MEDIUM_SIZE, handletextpad=0.1)
    plt.setp(legend.get_title(),fontsize=MEDIUM_SIZE)
    plt.savefig('./images/PDCD1-delfreq-revised.png', bbox_inches='tight', dpi=350)#linewidth=1

dotplt()

#141, 143 (143 not included before)

#VARIANT_DF SHOULD BE RENAMED PAM DF
#variant_df_path = './data/data/df/variant_df.csv'
#pam_df = pd.read_csv(variant_df_path)
#pam_df = pam_df[pam_df['genename'] == 'PDCD1'].dropna()
#maxcoord = max(pam_df['start'].tolist()) # because all are on reverse strand
#mincoord = min(pam_df['end'].tolist())


'''
coords = pam_df['start'].tolist().extend(pam_df['end'].tolist())
print(coords)

#var_df = pd.read_csv('data/variant/bed-and-tsv-before/PDCD1.tsv', sep='\t')
#print(var_df)

CDSpamsbed_path = './data/variant/CDSpams.bed'
CDSpamsbed = pd.read_csv(CDSpamsbed_path, sep='\t')
CDSpamsbed = CDSpamsbed[CDSpamsbed['#'] == 2]
CDSpamsbed =  CDSpamsbed[CDSpamsbed['start'] >= mincoord]
CDSpamsbed = CDSpamsbed[CDSpamsbed['end'] <= maxcoord]


genome = Fasta(genome_path)
seqs = []
for i in range(len(CDSpamsbed)):
    chrom = CDSpamsbed['chrom'].iloc[i]
    start = CDSpamsbed['start'].iloc[i] 
    end = CDSpamsbed['end'].iloc[i]
    if CDSpamsbed['strand'].iloc[i] == '+':
        seq = genome[chrom][start:end].seq
    else: # strand == '-'
        seq = genome[chrom][start:end].reverse.complement.seq
    seqs.append(seq)
'''