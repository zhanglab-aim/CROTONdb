'''
Create tabix file for database creation
050421 VL
'''

import os
import pysam
import logging
import pandas as pd
import argparse
from tqdm import tqdm

logger = logging.getLogger(__name__)
tabix_dir = './datavl/variant/tabix'
tsv_dir = './datavl/variant/tsv'

def make_tabixdirs():
    """
    Make /tabix/X, /tabix/Y, /tabix/1, /tabix/2 ... /tabix/22 directories in tabix_dir
        - So .gz and .gz.tbi files for CROTON_varpred.[chrom] dataframes can be stored by chromosome
    """
    chrom_lst = ['X'] #, 'Y'
    for i in range(22): chrom_lst.append(str(i+1))
    for chrom in chrom_lst:
        tabix_path = tabix_dir + '/%s' % (chrom)
        if not os.path.exists(tabix_path): # if path does not exist, create directory
            os.makedirs(tabix_path)
            print(tabix_path)

def combine_tsvs(chrom): # runs on subset of list of chromosomes by lstinx and num
    chrom = str(chrom)
    tsv_fp_lst = [os.path.join(path, name) for path, subdirs, files in os.walk(tsv_dir + '/%s'%(chrom)) for name in files]
    master_df = pd.DataFrame()
    for fp in tsv_fp_lst:
        df = pd.read_csv(fp, sep='\t')
        master_df = pd.concat([df, master_df]) ##chrom\tpos\tref\talt
    print('Successfully finished loop for chr%s!'%(chrom))
    master_df = master_df.rename(columns={'chrom':'#chrom', 'varid':'id'})
    master_df = master_df.sort_values(by=['#chrom', 'pos'])
    
    # chrom  pamid   start   end     strand  pos     id      ref     alt     AF_info ref_seq alt_seq ref_preds       alt_preds       abs_diffs       abs_diff
    # master_df = master_df[['chrom', 'pos', 'varid', 'ref', 'alt', 'pamid', 'genename', 'num', 'start', 'end', 
    #    'strand', 'inx', 'ref_seq', 'alt_seq', 'ref_preds', 'alt_preds', 'abs_diff', 'ref_del_frq', 'alt_del_frq',
    #    'diff_del_frq', 'ref_1ins', 'alt_1ins', 'diff_1ins', 'ref_1del', 'alt_1del', 'diff_1del', 'ref_onemod3', 'alt_onemod3',
    #    'diff_onemod3', 'ref_twomod3', 'alt_twomod3', 'diff_twomod3', 'ref_frameshift', 'alt_frameshift', 'diff_frameshift']]
    master_df.to_csv(tabix_dir + '/%s/CROTON_varpred_%s.tsv'%(chrom, chrom), sep='\t', index=False) # compression='gzip'
        
if __name__ == "__main__": # see sbatch_make_tabix.sh
    parser = argparse.ArgumentParser(description='combining tsvs in a vcf file')
    parser.add_argument('--chrom', help='chromosome to make vcf file for')
    args = parser.parse_args()
    make_tabixdirs()
    combine_tsvs(chrom=args.chrom)

#df_path = 'datavl/variant/tabix/22/CROTON_varpred_22.tsv'
#df = pd.read_csv(df_path, sep='\t')
#pri
# nt(df)

'''
# THEN RUN--> (1) gz compressions with bgzip, (2) make tabix files
(1) use bgzip to compress to .gz files
    for i in `seq 1 22`; do echo $i; bgzip -c datavl/variant/tabix/$i/CROTON_varpred_$i.tsv > datavl/variant/tabix/$i/CROTON_varpred_$i.gz; done
    for i in X Y; do echo $i; bgzip -c datavl/variant/tabix/$i/CROTON_varpred_$i.tsv > datavl/variant/tabix/$i/CROTON_varpred_$i.gz; done
(2) use tabix to make .gz.tbi files
    for i in `seq 1 22`; do echo $i; tabix -p vcf datavl/variant/tabix/$i/CROTON_varpred_$i.gz; done
    for i in X Y; do echo $i; tabix -p vcf datavl/variant/tabix/$i/CROTON_varpred_$i.gz; done
'''

def fetch_tabix_croton_feature_predictions(tabix_file_or_url, chrom=None, start=None, end=None): #https://github.com/aaronkw/humanbase-api/blob/4d5c6cb6eeefa902089d78c9acd435b2ffb498e5/humanbase/deepsea/utils.py#L253
    try:
        with pysam.TabixFile(tabix_file_or_url, encoding='utf-8') as t:
            if chrom:
                tabix_result = list(
                    t.fetch(
                        'chr{}'.format(chrom),
                        start,
                        end,
                        parser=pysam.asTuple(),
                    )
                )
            else:
                tabix_result = list(t.fetch(parser=pysam.asTuple()))
    except ValueError:
        logger.debug('No data found for positions {}-{}'.format(start, end))
        return []

    if not len(tabix_result):
        logger.debug('No data found for positions {}-{}'.format(start, end))
        return []

    if len(t.header):
        headers = t.header[0].split('\t')
        return pd.DataFrame(tabix_result, columns=headers)

    return pd.DataFrame(tabix_result)

#bgzip -c datavl/variant/tabix/22/CROTON_varpred_22.tsv > datavl/variant/tabix/22/CROTON_varpred_22.vcf.gz
#tabix -p vcf datavl/variant/tabix/22/CROTON_varpred_22.vcf.gz

# should be chr then position for vcf format
# use bgzip
# put a '#' in the header
# when adding to server, make it look as close to previous file versions as possible, add one col at a time
# create .bash file going through all scripts by chrom
