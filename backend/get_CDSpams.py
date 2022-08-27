'''
Get all PAM sites in the coding region (by chrom) then intersect with .vcf files describing SNVs
VRL
    - Created May 2021? 
    - Edited 072221
'''

import gffutils
import pandas as pd
import numpy as np
import json
from pyfaidx import Fasta
import matplotlib.pyplot as plt
import seaborn as sns
import bisect
import intervaltree
from tqdm import tqdm

pam_df_path = './datavl/df/pam_df.csv'
genome_path = './data/data/genome/GRCh38.primary_assembly.genome.fa' #/mnt/ceph/users/zzhang/human_SNP/Gnomad_V2_Exomes_hg38_liftover/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz
CDSpamsbed_path = './datavl/variant/CDSpams.bed'
db_path = '/mnt/ceph/users/zzhang/decan/data/genomes/gencode/gencode.v35.annotation.gff3.gz.decan_db'
db = gffutils.FeatureDB(db_path)
gn_dict = {g.attributes['gene_name'][0]: g.id for g in db.features_of_type('gene')}

CDSpamsbyChrom_dir = './datavl/variant/CDSpams-byChrom/'

#######################
### Make pam_df.csv ###
#######################
# Contains columns: 'genename', 'chrom', 'start', 'end', 'strand', 'seq', 'pams', 'rc_pams'

def get_CDS_df():
    """
    Yields
    ------
    Dataframe with 'genename', 'chrom', 'start', 'end, 'strand' columns for all CDS detected in db
    - Save at pam_df_path defined at top of the file
    """
    global db
    CDS_feature_lst = list(db.all_features(featuretype='CDS')) #774993
    genenames, chroms, starts, ends, strands = [], [], [], [], []
    for inx in range(len(CDS_feature_lst)):
        feature = CDS_feature_lst[inx].astuple()
        genename = json.loads(feature[9])['gene_name'][0]
        chrom = feature[1]
        start = feature[4]
        end = feature[5]
        strand = feature[7]
        genenames.append(genename)
        chroms.append(chrom)
        starts.append(start)
        ends.append(end)
        strands.append(strand)

    features_dict = {'genename':genenames, 'chrom':chroms, 'start': starts, 'end': ends, 'strand': strands}

    df = pd.DataFrame(features_dict)
    df.to_csv(pam_df_path, index=False)

def get_seq(): 
    """
    Add 'seq' column to pam_df by using 'start', 'end', 'chrom' columns + pyfaidx
    """
    df = pd.read_csv(pam_df_path)
    
    # Only keep chr1-22, chrX, CDS in pam_df
    #chrom_lst = ['chrX', 'chrY'] 
    chrom_lst = ['chrX']
    for i in range(22): chrom_lst.append('chr%i'%(i+1))
    df = df.loc[df['chrom'].isin(chrom_lst)] #774980

    # Get sequences and add to dataframe
    genome = Fasta(genome_path)
    seqs = []
    for i in range(len(df)):
        chrom = df['chrom'].iloc[i]
        start = df['start'].iloc[i] - 1
        end = df['end'].iloc[i]
        if df['strand'].iloc[i] == '+':
            seq = genome[chrom][start:end].seq
        else: # strand == '-'
            seq = genome[chrom][start:end].reverse.complement.seq
        seqs.append(seq)
    df['seq'] = seqs
    df.to_csv(pam_df_path, index=False)

# add PAM columns to pam_df.csv

def findall(p, s):
    '''
    Yields
    ------
    All positions of p in the string s
    '''
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i + 1)

def reverse_complement(seq): # get reverse complement of sequence
    letter_match = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    reverse_complement = "".join(letter_match[b] for b in reversed(seq))
    return reverse_complement

def find_pam_coords(seq, start, end):
    '''
    Get pam coords of pams for both + and - versions of a sequence 

    Parameters
    ----------
    seq: string 
        DNA sequence
    start: int
        Start coordinate of sequence
    end: int
        End coordinate of sequence
    
    Yields
    ------
    (List of coords of pams on 'seq', List of coords of pams on reverse complement of 'seq')
        - Pam coordinates correspond with the 'N' of 'NGG'
    
    Examples
    --------
    >>> print(find_pam_coords('ACCGCCAGGG', 0, 9))
        Output: ([6, 7], [3, 6]) -->
            seq:                              rc_seq:
            A  C  C  G  C  C  A  G  G  G  //  C  C  C  T  G  G  C  G  G  T
            0  1  2  3  4  5 *6 *7  8  9  //  0  1  2 *3  4  5 *6  7  8  9
    '''

    possible_pams = ['AGG', 'CGG', 'GGG', 'TGG'] # 'NGG' PAMs
    rc_seq = reverse_complement(seq)
    pam_inxs = [] 
    rc_pam_inxs = [] # reverse complement
    for i in range(len(possible_pams)):
        pam = possible_pams[i]
        pam_inx_lst = [i for i in findall(pam, seq)]
        pam_inxs.extend(pam_inx_lst)
        rc_pam_inx_lst = [i for i in findall(pam, rc_seq)]
        rc_pam_inxs.extend(rc_pam_inx_lst)
    
    coords = [x + start for x in pam_inxs] #assuming start inclusive
    rc_coords = [end - x for x in rc_pam_inxs] #and end exclusive (why -1)

    return coords, rc_coords

def get_pam_coords():
    """
    Add columns ['pams', 'rc_pams']  to pam_df by using find_pam_coords function [rc = rverse complement]
        - Columns contain lists (in str form) of indices of pams ('N' of 'NGG')
    
    Yields 
    ------
        Final pam_df with columns ['genename', 'chrom', 'start', 'end', 'strand', 'seq', 'pams', 'rc_pams']
    
    Example
    -------
    >>> First three rows --> 
         genename chrom start   end strand                       seq                     pams                  rc_pams
            OR4F5  chr1 65565 65573      +                  TGAAGAAG                       []                       []
            OR4F5  chr1 69037 70008      +  TAACTGCAGAGGCTATTTCCT...  [69046, 69130, 69176...  [69909, 69881, 69847...
            OR4F5  chr1 69091 70008      +  TGGTGACTGAATTCATTTTTC...  [69130, 69176, 69352...  [69909, 69881, 69847...
    """
    df = pd.read_csv(pam_df_path)
    df['pams'] = np.nan
    df['rc_pams'] = np.nan
    df['pams'] = df['pams'].astype('object') # so list can be saved to df
    df['rc_pams'] = df['rc_pams'].astype('object')
    for i in range(len(df)):
        seq = df['seq'].iloc[i]
        start = df['start'].iloc[i] 
        end = df['end'].iloc[i]
        strand = df['strand'].iloc[i]
        if start - end != 0: #nan 1266290 1266290
            if strand == '+': coords, rc_coords = find_pam_coords(seq, start, end)
            else: coords, rc_coords = find_pam_coords(reverse_complement(seq), start, end) # WAS A TYPO HERE coord instead of coords
            df.at[i, 'pams'] = coords 
            df.at[i, 'rc_pams'] = rc_coords
    
    df.to_csv(pam_df_path, index=False)

########################
### Make CDSpams.bed ###
########################
#pd.set_option("max_rows", None)

def make_CDSpamsbed():
    """
    Take pam_df, drop nan for rows that start and end with nan; drop rows w/o any pams
    
    Notes [Procedure]
    -----
    For each gene: 
        - Combine lists in 'pams' column into one list; combine lists in 'rc_pams' column into one list [keep unique]
        - Sum of lenghts of these two is the total number of pams for a gene
        - Get start and end coordinates from pam coord
            - Sort in ascending order by 'start'
              --> Cut sites should be in order reading left to right on + strand (if 2 pam sites same 60 bp sequence, reverse complement first b/c pam on left of cutsite)
        - Then add columns: '#': chr, 'genename', 'num': numerical identification for PAM within gene, 'pamid': unique id for pam, genename|num
    
    Returns
    ------
    The file CDSpam.bed at CDSpamsbed_path defined at top of this script
    
    Example 
    -------
    >>> First couple of rows of CDSpam.bed
        start       end strand  # genename  num     pamid
        69013     69073      +  1    OR4F5    1   OR4F5|1
        69024     69084      +  1    OR4F5    2   OR4F5|2
        69031     69091      -  1    OR4F5    3   OR4F5|3
        69058     69118      +  1    OR4F5    4   OR4F5|4
        69079     69139      +  1    OR4F5    5   OR4F5|5
    """
    
    pam_df = pd.read_csv(pam_df_path).dropna() # 394 all start & end with nan
    pam_df = pam_df[(pam_df['pams'] != '[]') | (pam_df['rc_pams'] != '[]')].reset_index(drop=True) #12033 rows
    pam_df = pam_df[pam_df['chrom']!='chrY']
    CDSpams_chr, CDSpams_gn, CDSpams_num = [], [], []
    CDSpamsbed_df = pd.DataFrame()
    #genome = Fasta(genome_path)
    
    skipped_genes = []
    for genename in tqdm(list(pam_df['genename'].unique())):
        gn_start, gn_end, gn_strand, gn_pam = [], [], [], [] 
        #gn_seq = []
        pam_df_gn = pam_df[pam_df['genename'] == genename]
        if not len(pam_df_gn['chrom'].unique()) == 1:
            skipped_genes.append(genename)
            continue
        pam_lst = list(set([a for b in pam_df_gn['pams'].str.strip('][').str.split(', ') for a in b])) # different rows for the same gene have overlapping ranges
        rc_pam_lst = list(set([a for b in pam_df_gn['rc_pams'].str.strip('][').str.split(', ') for a in b]))
        if '' in pam_lst: pam_lst.remove('')
        if '' in rc_pam_lst: rc_pam_lst.remove('')
        
        num_pams = len(pam_lst) + len(rc_pam_lst)
        chrom = pam_df_gn['chrom'].iloc[0]
        chr_lst = num_pams * [chrom]
        gn_lst = num_pams * [genename]
        num_lst = [i + 1 for i in range(num_pams)]
        CDSpams_chr += chr_lst
        CDSpams_gn += gn_lst
        CDSpams_num += num_lst

        for pam in rc_pam_lst: # '-' cases have a lower ID number because PAM on the left side of the cutsite
            pam = int(pam)
            start = pam - 27
            end = pam + 33
            gn_start.append(start)
            gn_end.append(end)
            gn_strand.append('-')
            #seq = genome[chrom][start:end].reverse.complement.seq
            #gn_seq.append(seq)
        
        for pam in pam_lst:
            pam = int(pam)
            start = pam - 33
            end = pam + 27
            gn_start.append(start)
            gn_end.append(end)
            gn_strand.append('+')
            #seq = genome[chrom][start:end].seq
            #gn_seq.append(seq)
        
        gene_dict = {'start':gn_start, 'end':gn_end, 'strand':gn_strand} # , 'pam':pam_lst + rc_pam_lst, 'seq':gn_seq
        gene_df = pd.DataFrame(gene_dict).sort_values(by='start')
        CDSpamsbed_df = CDSpamsbed_df.append(gene_df)

    CDSpamsbed_df['#'] = CDSpams_chr
    CDSpamsbed_df['#'] = CDSpamsbed_df['#'].str.replace('chr', '')
    CDSpamsbed_df['genename'] = CDSpams_gn
    CDSpamsbed_df['num'] = CDSpams_num
    CDSpamsbed_df['pamid'] = CDSpamsbed_df['genename'] + '|' + CDSpamsbed_df['num'].astype(str)
    
    #print(CDSpamsbed_df.shape)
    print(skipped_genes)
    CDSpamsbed_df.to_csv(CDSpamsbed_path, index=False, sep='\t') #5368021 rows x 7 columns

#make_CDSpamsbed()

##################################
## Make X.bed, 1.bed, 2.bed ... ##
##################################

def subset_CDSpams(): 
    '''
    Make X.bed, 1.bed, 2.bed ... in CDSpamsbyChrom_dir by taking subsets of CDSpams.bed
    '''
    CDSpamsbed_df = pd.read_csv(CDSpamsbed_path, sep='\t')
    chrom_lst = ['X'] # ['X', 'Y']
    for i in range(22): chrom_lst.append(str(i+1))
    for chrom in chrom_lst:
        print(chrom)
        chrom_df = CDSpamsbed_df[CDSpamsbed_df['#'].astype(str) == chrom] # before some of chrom22 has integer 22 as a col value
        chrom_df = chrom_df[['#', 'start', 'end', 'strand', 'pamid', 'genename', 'num']]
        chrom_df['#'] = 'chr' + chrom_df['#'].astype(str) # add 'chr' to chromosome '#' column
        chrom_df.to_csv(CDSpamsbyChrom_dir + '%s.bed' % (str(chrom)), sep='\t', index=False)

#subset_CDSpams()

if __name__ == "__main__":
    # Make pam_df.csv
    print('Get df')
    get_CDS_df()
    print('Get seq')
    get_seq()
    print('Get pam coords')
    get_pam_coords()

    # Make CDSpams.bed
    print('Make CDSpamsbed')
    make_CDSpamsbed()
    print('Subset CDSpams')
    subset_CDSpams()

################################################################
## Make X.intersect.bed, 1.intersect.bed, 2.intersect.bed ... ##
################################################################
# - Decided to not do the Y chromosome because of repeated genes between X and Y

'''
$ cd datavl/variant 
'''

###### original version (no population data) ######

# $ for i in `seq 1 22`; do echo $i; bedtools intersect -a CDSpams-byChrom/$i.bed -b /mnt/home/zzhang/ceph/human_SNP/hg38_All_20180418.vcf-byChrom/$i.vcf -wo > CDSpams-byChrom/$i.intersect.bed; done
# $ bedtools intersect -a CDSpams-byChrom/X.bed -b /mnt/home/zzhang/ceph/human_SNP/hg38_All_20180418.vcf-byChrom/X.vcf -wo > CDSpams-byChrom/X.intersect.bed

######*** for the 1000 genome files (with population data) ######***
''' To create .vcf files to intersect .bed files with (in 'human_SNP' NOT 'croton' folder) 

--> NOTE: 09/10/21 Frank filtered out variants w/o varid (there's prob a different way to make these now) smth like

$ for i in `seq 1 22`; do echo $i; tabix /mnt/ceph/users/zzhang/human_SNP/Gnomad_V2_Exomes_hg38_liftover/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz chr$i > /mnt/ceph/users/vli/human_SNP/Gnomad_V2_Exomes_hg38_liftover/chr$i.vcf; done 
$ for i in X Y; do echo $i; tabix /mnt/ceph/users/zzhang/human_SNP/Gnomad_V2_Exomes_hg38_liftover/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz chr$i > /mnt/ceph/users/vli/human_SNP/Gnomad_V2_Exomes_hg38_liftover/chr$i.vcf; done 
$ for i in `seq 1 22`; do echo $i; echo "##fileformat=VCFv4.2" | cat - /mnt/ceph/users/vli/human_SNP/Gnomad_V2_Exomes_hg38_liftover/chr$i.vcf > /mnt/ceph/users/vli/human_SNP/Gnomad_V2_Exomes_hg38_liftover2/chr$i.vcf; done 
$ for i in X Y; do echo $i; echo "##fileformat=VCFv4.2" | cat - /mnt/ceph/users/vli/human_SNP/Gnomad_V2_Exomes_hg38_liftover/chr$i.vcf > /mnt/ceph/users/vli/human_SNP/Gnomad_V2_Exomes_hg38_liftover2/chr$i.vcf; done 
'''
# EDITED 09/10/21
# $ for i in `seq 1 22`; do echo $i; bedtools intersect -a CDSpams-byChrom/$i.bed -b /mnt/home/zzhang/ceph/human_SNP/Gnomad_V2_Exomes_hg38_liftover/by_chrom/chr$i.vcf -loj > CDSpams-byChrom/$i.intersect.bed; done 
# $ bedtools intersect -a CDSpams-byChrom/X.bed -b /mnt/home/zzhang/ceph/human_SNP/Gnomad_V2_Exomes_hg38_liftover/by_chrom/chrX.vcf -loj > CDSpams-byChrom/X.intersect.bed

#/mnt/home/zzhang/ceph/human_SNP/Gnomad_V2_Exomes_hg38_liftover/by_chrom for filtered vcfs

###########################################
## function to get CDS for specific gene ##
###########################################

def get_CDS_for_gene(genename):
    """
    
    """
    global db, gn_dict
    gid = gn_dict[genename]
    gene = db[gid]
    order_by = 'start' if gene.strand == '+' else 'end'
    cdss = list(db.children(gid, featuretype='CDS', order_by=order_by))
    it = intervaltree.IntervalTree()
    for cds in cdss:
        try:
            it[cds.start:cds.end] = 1
        except ValueError:
            pass
    #cds_array = sorted(list(set([cds.start if gene.strand=='+' else cds.end for cds in cdss]))) # for lookup
    it.merge_overlaps()
    cds_array = sorted([x.begin if gene.strand=='+' else x.end for x in it])
    return cds_array

#pam_df = pd.read_csv(pam_df_path, nrows=100).dropna()
#print(pam_df)