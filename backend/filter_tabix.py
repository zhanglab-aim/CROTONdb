'''
Find CROTON predicted high impact variants, with 
(a) consistent high allele frequency (ex. >5%), and 
(b) highly variable allele frequency among genetic ancestries (ex. some >5%, some <1)?
VL 092521
'''
import os
import ast
import pysam
import argparse
import pandas as pd

# write down where files are stored in web dev google doc
# tune parameters , look into it
# 'blood', 'T lymphocyte' https://hb.flatironinstitute.org/module/overview/ human base couple of blood networks, lymphocyte
# annotate this more

AF_lst = ['AF_afr', 'AF_amr', 'AF_asj', 'AF_eas', 'AF_fin', 'AF_nfe', 'AF_oth', 'AF_sas'] # not sex specific
out_dir = './datavl/variant/filt-fs-AF'

def make_filtdirs():
    """
    Make directories to store the files created in this script

    Yield Directories (if they do not exist)
    ---------------------------------------- 
    >>> './datavl/variant/filt-fs-AF/filt-fs' : get_high_impact_df()
        - files with a high abs diff btw ref and alt frameshift (fs) freq
    >>> './datavl/variant/filt-fs-AF/consis_highAF' : consis_highAF()
        - files with high fs and consistently high AF
    >>> './datavl/variant/filt-fs-AF/variable_AF' : variable_AF()
        - files with high fs and variable AF
    """ 
    if not os.path.exists(out_dir + '/filt-fs'):
        os.makedirs(out_dir + '/filt-fs')
        print(out_dir + '/filt-fs')
    if not os.path.exists(out_dir + '/consis-highAF'): 
        os.makedirs(out_dir + '/consis-highAF')
        print(out_dir + '/consis-highAF')
    if not os.path.exists(out_dir + '/variable-AF'): 
        os.makedirs(out_dir + '/variable-AF')
        print(out_dir + '/variable-AF')

def get_highimpact_df(chrom, start=None, end=None, diff_annot='frameshift', diff_filter=0.2): 
    """
    Get dataframes of variants with absolute predicted frequencies > diff_filter
        Tabix files from directory '.datavl/variant/tabix'
    
    Parameters 
    ----------
    chrom: int or str
    start: int or None
    end: int or None
        chrom, start, and end define the region to fetch from Tabix files
        when start = None and end = None entire Tabix file will be fetched
    diff_col : str
        the annotation/task out of 6 predicted CROTON values to compute difference on
    diff_filter: float
        lower limit of abs fs freq differences = abs(ref - alt)
    
    Yield
    -----
    Dataframes with absolute frameshift differences above fs_filter 
        - Note: rows without variants have been filtered out as well (varid = '.')
    """
    assert diff_annot in ('frameshift', '1ins', '1del', 'onemod3', 'twomod3')
    tabixfile = pysam.TabixFile('datavl/variant/tabix/%s/CROTON_varpred_%s.gz' % (str(chrom), str(chrom)))
    col_lst = tabixfile.header[0].split('\t')
    df = pd.DataFrame(columns=col_lst)
    tabix_result = list(tabixfile.fetch('chr'+str(chrom), start, end))
    print(len(tabix_result)) # number of rows in fetched tabix file region
    
    col_dict = {i:x for x,i in enumerate(col_lst)}
    ref_col = col_dict.get(f"ref_{diff_annot}")
    alt_col = col_dict.get(f"alt_{diff_annot}")

    for i in range(len(tabix_result)):
        result = tabix_result[i].split('\t')
        if (float(result[-1]) != -1) and (abs(float(result[ref_col]) - float(result[alt_col])) > diff_filter): 
            # Note: last columns/entries in list are 'ref_frameshift' (-2) and 'alt_frameshift' (-1), respectively
            # filter for high abs fs freq diff, remove rows with no var (alt_fs == '-1')
            df = df.append(pd.Series(result, index=col_lst), ignore_index=True)
    
    print('Finished filtering df!')
    print(len(df))
    df.to_csv((out_dir + '/filt-fs/%s.tsv' % str(chrom)), index=False, sep='\t')
    return df

def consis_highAF(chrom, min=0.05, start=None, end=None, diff_annot='frameshift', diff_filter=0.2): # (a)
    """
    Get dataframes with high fs freq and consistently high allele frequency (AF)
    
    Parameters 
    ----------
    chrom, start, end, fs_filter --> see get_highimpact_df() 
    min: float
        Lower limit of what is considered a high AF (noninclusive)
    
    Yield
    -----
    Dataframes with entries where ALL AFs in the AF_list are > min
    """

    try: df = pd.read_csv(out_dir + '/filt-fs/%s.tsv' % str(chrom), sep='\t')
    except: df = get_highimpact_df(chrom, start, end, diff_annot=diff_annot, diff_filter=diff_filter)
    for i in range(len(df)):
        AFdict = ast.literal_eval(df['AF_info'].loc[i])
        filt_AFdict = {key: value for key, value in AFdict.items() if key in AF_lst and float(value) > min} 
        # dictionary contains all AFs > min
        if len(filt_AFdict.keys()) != 8: # if not all 8 AFs > min, drop row
            df = df.drop(labels=i, axis=0)
    
    print('Saving consistent highAF...')
    print(len(df))
    df.to_csv((out_dir + '/consis-highAF/%s.tsv' % str(chrom)), index=False, sep='\t')
    
    return df
    
def variable_AF(chrom, lower=0.01, upper=0.05, start=None, end=None, diff_annot='frameshift', diff_filter=0.2): #b
    """
    Get dataframes with high fs freq and variable allele frequency (AF)
    
    Parameters 
    ----------
    chrom, start, end, fs_filter --> see get_highimpact_df() 
    lower: float
        Upper limit of what is considered a 'low' AF
    upper: float
        Lower limit of what is considered a 'high' AF
    
    Yield
    -----
    Dataframes with entries where some AFs are below the 'lower'limit 
        and some are above the 'upper'limit
    """
    try: df = pd.read_csv(out_dir + '/filt-fs/%s.tsv' % str(chrom), sep='\t')
    except: df = get_highimpact_df(chrom, start, end, diff_annot=diff_annot, diff_filter=diff_filter)
    for i in range(len(df)):
        AFdict = ast.literal_eval(df['AF_info'].loc[i])
        filt_AFdict = {key: value for key, value in AFdict.items() if key in AF_lst} # dict only contains keys of AFs in AF_lst
        below_lower = {key: value for key, value in filt_AFdict.items() if float(value) < lower} # dict with AF values < lower
        above_upper = {key: value for key, value in filt_AFdict.items() if float(value) > upper} # dict with AF values > upper
        if (len(below_lower.keys()) == 0) or (len(above_upper.keys()) == 0): 
            # if row does not have AFs below 'lower' or above 'upper',)drop it
            df = df.drop(labels=i, axis=0)
        
    print('Saving variable AF...')
    print(len(df))
    df.to_csv((out_dir + '/variable-AF/%s.tsv' % str(chrom)), index=False, sep='\t')
    
    return df

'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Making filtered fs and AF dataframes')
    parser.add_argument('--chrom', help='Chromosome with which to make filtered dataframes')
    args = parser.parse_args()
    make_filtdirs()
    print('Making files for chr', args.chrom)
    consis_highAF(chrom=args.chrom)
    variable_AF(chrom=args.chrom)

#for i in `seq 1 22`; python src/variant/analyze_tabix.py --chrom $i; done
#python src/variant/analyze_tabix.py --chrom X
'''

def combine_files(dirname):
    """
    Combine the files in a directory (in this case tsvs from consis_highAF and variable_AF)

    Parameters 
    ----------
    dirname: str
        Name of directory containing dataframes to combine
    
    Yield
    -----
    Merged dataframe that combines all dataframes in dirname directory
    """
    
    file_list = os.listdir(dirname)
    list_of_dataframes = []
    for filename in file_list:
        list_of_dataframes.append(pd.read_csv((dirname + filename), sep='\t'))
    
    merged_df = pd.concat(list_of_dataframes)
    merged_df.to_csv((dirname + 'merged.tsv'), index=False, sep='\t')
    print(merged_df)
    
    return merged_df

make_filtdirs()
for chrom in range(1, 23):
    chrom = str(chrom)
    print('Making files for chr', chrom)
    consis_highAF(chrom=chrom, diff_annot="1ins", diff_filter=0.2)
    variable_AF(chrom=chrom, diff_annot="1ins", diff_filter=0.2)
combine_files('./datavl/variant/filt-fs-AF/consis-highAF/')
combine_files('./datavl/variant/filt-fs-AF/variable-AF/')

'''
FOR PARAMETERS: fs_filter = 0.2, min = 0.05lower=0.01, upper=0.05
    ran make_filedirs(), consis_highAF, variable_AF for all chroms
    When combine: variable-AF 576 entries, consis-highAF 363 entries

    Making files for chr 1
    889228
    Finished filtering df!
    976
    Saving consistent highAF...
    22
    Saving variable AF...
    83

    Making files for chr 2
    542557
    Finished filtering df!
    506
    Saving consistent highAF...
    22
    Saving variable AF...
    42

    Making files for chr 3
    441861
    Finished filtering df!
    460
    Saving consistent highAF...
    17
    Saving variable AF...
    23

    Making files for chr 4
    275948
    Finished filtering df!
    246
    Saving consistent highAF...
    12
    Saving variable AF...
    9

    Making files for chr 5
    344308
    Finished filtering df!
    367
    Saving consistent highAF...
    15
    Saving variable AF...
    10

    Making files for chr 6
    402390
    Finished filtering df!
    483
    Saving consistent highAF...
    23
    Saving variable AF...
    30

    Making files for chr 7
    413354
    Finished filtering df!
    505
    Saving consistent highAF...
    21
    Saving variable AF...
    29

    Making files for chr 8
    299028
    Finished filtering df!
    323
    Saving consistent highAF...
    14
    Saving variable AF...
    18

    Making files for chr 9
    366244
    Finished filtering df!
    364
    Saving consistent highAF...
    20
    Saving variable AF...
    19

    Making files for chr 10
    316856
    Finished filtering df!
    336
    Saving consistent highAF...
    15
    Saving variable AF...
    30

    Making files for chr 11
    548121
    Finished filtering df!
    592
    Saving consistent highAF...
    28
    Saving variable AF...
    35

    Making files for chr 12
    406529
    Finished filtering df!
    460
    Saving consistent highAF...
    30
    Saving variable AF...
    20

    Making files for chr 13
    121376
    Finished filtering df!
    126
    Saving consistent highAF...
    1
    Saving variable AF...
    2

    Making files for chr 14
    287089
    Finished filtering df!
    272
    Saving consistent highAF...
    8
    Saving variable AF...
    9

    Making files for chr 15
    293042
    Finished filtering df!
    308
    Saving consistent highAF...
    11
    Saving variable AF...
    32

    Making files for chr 16
    491361
    Finished filtering df!
    631
    Saving consistent highAF...
    16
    Saving variable AF...
    34

    Making files for chr 17
    564652
    Finished filtering df!
    772
    Saving consistent highAF...
    10
    Saving variable AF...
    44

    Making files for chr 18
    115170
    Finished filtering df!
    132
    Saving consistent highAF...
    9
    Saving variable AF...
    2

    Making files for chr 19
    751753
    Finished filtering df!
    1035
    Saving consistent highAF...
    40
    Saving variable AF...
    43

    Making files for chr 20
    236531
    Finished filtering df!
    255
    Saving consistent highAF...
    10
    Saving variable AF...
    16

    Making files for chr 21
    98643
    Finished filtering df!
    164
    Saving consistent highAF...
    2
    Saving variable AF...
    10

    Making files for chr 22
    241706
    Finished filtering df!
    305
    Saving consistent highAF...
    10
    Saving variable AF...
    19

    Making files for chr X
    271239
    Finished filtering df!
    251
    Saving consistent highAF...
    7
    Saving variable AF...
    17
'''
