'''
Make bed files, one for each gene

Author: VRL

This step is not necessary; only kept here for legacy reasons and to continue work from next pipeline step (make_tsvs). 
Adds ~120min to finish when running sequentially. FZZ 20220827
'''

import pandas as pd
import argparse
import os
import numpy as np
from tqdm import tqdm
from . import configs


################
### Get beds ###
################
#gen or ccb

def make_BEDdirs():
    """
    Make /bed/X, /bed/Y, /bed/1, /bed/2 ... /bed/22 directories in bed_dir
        - So .bed files for genes can be stored according to their chromosome
    """
    chrom_lst = ['X'] # 'Y'
    for i in range(22): chrom_lst.append(str(i+1))
    for chrom in chrom_lst:
        bedchr_path = configs.GENE_BED_DIR + '/%s' % (chrom)
        if not os.path.exists(bedchr_path): # if path does not exist, create directory
            os.makedirs(bedchr_path)
            #print(bedchr_path)

def get_GENE_beds(chrom):
    """
    Create bed files for all genes in a chromosome

    Returns
    -------
    list : a list of genes without any variants
    """
    chrom = str(chrom)
    print('Making bed files for chr' + chrom)
    novar_name = []
    bed_df = pd.read_csv(configs.CDS_PAM_DIR +'/%s.intersect.bed' % (chrom), sep='\t', 
        names=['chrom', 'start', 'end', 'strand', 'pamid', 'genename', 'num', 'chrom_', 'pos', 'varid', 'ref', 'alt', '.', '..', 'AF_info']) 
    print('Loaded intersect.bed!')
    genenames_chrom = sorted(bed_df['genename'].unique().tolist())
    for genename in tqdm(genenames_chrom):
        #print(genename)
        gn_df = bed_df[bed_df['genename'] == genename].reset_index(drop=True)
        for i in range(len(gn_df)): 
            row_info = gn_df['AF_info'].iloc[i]
            row_info = row_info.split(';')
            row_info = [entry for entry in row_info if '=' in entry] # Ex: some columns are only 'segdup'
            row_dict = dict(map(lambda x: x.split('='), row_info))
            row_dict_AF = {k: row_dict[k] for k in configs.GNOMAD_AF_ENTRIES if k in row_dict}
            gn_df['AF_info'][i] = str(row_dict_AF)
    
        if len(gn_df.index) != 0: gn_df.to_csv(configs.GENE_BED_DIR+ '/%s/%s.bed' % (chrom, genename), index=False, sep='\t')
        else: novar_name.append(genename)
    #print(novar_name)
    return novar_name

'''
if __name__ == "__main__": # see sbatch_bed.sh
    parser = argparse.ArgumentParser(description='Making beds')
    parser.add_argument('--chrom', help='Chromosome with which to make beds')
    args = parser.parse_args()
    make_beddirs()
    get_beds(chrom=args.chrom)
'''
def get_beds_before(lstinx, num=8):
    """
    Get all genenames with pams from pam_df, sort genenames, and then make [genename].bed files
    
    Parameters [so can split up work on multiple slurm jobs]
    ----------
    num: int
        Number of 'buckets' to split sorted genenames list into [will have a list of lists]
    lstinx : int in the interval [0, num)
        Index of 'bucket' from split genenames list to run on [inx of a list in list of lists]
    
    Yield
    -----
    ([genename].bed) files in (bed_dir + /[chromosome]) directories
    
    Example
    -------
    >>> first couple of rows of PDCD1.bed
        chrom      start        end strand      pamid genename  num  chrom_        pos         varid        ref  alt  . ..                                               info  1
            2  241851028  241851088      +    PDCD1|1    PDCD1    1       2  241851031   rs376045669          C  A,T  .  .  RS=376045669;RSPOS=241851031;dbSNPBuildID=138;...  1
            2  241851028  241851088      +    PDCD1|1    PDCD1    1       2  241851040  rs1383867797          G    A  .  .  RS=1383867797;RSPOS=241851040;dbSNPBuildID=151...  1
            2  241851028  241851088      +    PDCD1|1    PDCD1    1       2  241851050  rs1382434393          G    A  .  .  RS=1382434393;RSPOS=241851050;dbSNPBuildID=151...  1
            2  241851028  241851088      +    PDCD1|1    PDCD1    1       2  241851055   rs770917505          C    T  .  .  RS=770917505;RSPOS=241851055;dbSNPBuildID=144;...  1
    """
    pam_df = pd.read_csv(pam_df_path) 
    pam_df = pam_df.dropna() # 394 all start = end with nan
    print(pam_df)
    pam_df = pam_df[(pam_df['pams'] != '[]') | (pam_df['rc_pams'] != '[]')].reset_index(drop=True)
    genenames = sorted(list(pam_df['genename'].unique()))
    genenames = [genenames[i::num] for i in range(num)][lstinx] #genenames = ['ACE2', 'CTLA4', 'CCR5', 'PDCD1']
    novar_name = []
    for genename in genenames:
        print(genename)
        chrom = pam_df[pam_df['genename'] == genename]['chrom'].tolist()[0]
        bed_df = pd.read_csv('./datavl/variant/CDSpams-byChrom/%s.intersect.bed' % (chrom.replace('chr', '')), 
            sep='\t', names=['chrom', 'start', 'end', 'strand', 'pamid', 'genename', 'num', 'chrom_', 'pos', 'varid', 'ref', 'alt', '.', '..', 'info', '1'])
        bed_df = bed_df[bed_df['genename'] == genename]
        if len(bed_df.index) != 0: bed_df.to_csv(bed_dir + '/%s/%s.bed' % (chrom.replace('chr', ''), genename), index=False, sep='\t')
        else: novar_name.append(genename)
    print(novar_name)
