'''
Run CROTON make tsvs

'''

import pandas as pd
from pyfaidx import Fasta
import tensorflow as tf
from tensorflow.python.keras.models import load_model
import numpy as np
import os
import argparse
import configs

#tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

model_fp = "./outputs/forecast_freqs-4/train.bak6.wsf6/bestmodel.h5"
genome_path = './data/data/genome/GRCh38.primary_assembly.genome.fa'
bed_dir = './datavl/variant/bed'
tsv_dir = './datavl/variant/tsv'

#pd.set_option('display.max_rows', 500)

def make_tsvdirs():
    """
    Make /tsv/X, /tsv/Y, /tsv/1, /tsv/2 ... /tsv/22
        - So .tsv files for genes can be stored according to their chromosome
    """
    chrom_lst = ['X', 'Y']
    for i in range(22): chrom_lst.append(str(i+1))
    for chrom in chrom_lst:
        tsvchr_path = tsv_dir + '/%s' % (chrom)
        if not os.path.exists(tsvchr_path): # if path does not exist, create directory
            os.makedirs(tsvchr_path)
            print(tsvchr_path)

def one_hot_encode(seq, base_map): 
    '''
    One hot encode DNA sequences
    
    Parameters 
    ----------
    seq: str
        sequence to be one hot encoded
    base_map: str where len(bas_map) = 4
        base_map[inx of letter] --> array converted into
        base_map[0] --> [1, 0, 0, 0]
        base_map[1] --> [0, 1, 0, 0]
        base_map[2] --> [0, 0, 1, 0]
        base_map[3] --> [0, 0, 0, 1]
    
    Returns
    -------
    One hot encoded 4 by n array where n = len(seq)

    '''
    seq = seq.upper()
    mapping = dict(zip(base_map, range(4)))
    split_seq = seq.split('N')
    map_seq = [mapping[i] for i in split_seq[0]]
    final_seq = np.eye(4)[map_seq]
    for n in range(len(split_seq) - 1):
        N_arr = np.array([0.25, 0.25, 0.25, 0.25])
        final_seq = np.vstack((final_seq, N_arr))
        map_seq = [mapping[i] for i in split_seq[n + 1]]
        final_seq_ = np.eye(4)[map_seq]
        final_seq = np.vstack((final_seq, final_seq_))
    return final_seq

def reverse_complement(seq):
    '''Get reverse complement of DNA sequence'''
    letter_match = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    reverse_complement = "".join(letter_match[b] for b in reversed(seq))
    return reverse_complement

def get_varscore(bed_fp, model):
    """
    Run CROTON on 60 bp sequences in bed_df at bed_fp (in bed_dir)
    
    Parameters
    ----------
    bed_fp: filepath to a bed_df [made in make_beds.py]
        - ex:  './data/variant/bed/2/PDCD1.bed'
    model: a Keras model instance
        - i.e. output of load_model function
    """
    df = pd.read_csv(bed_fp, sep='\t')
    df = df[['chrom', 'pos', 'start', 'end', 'strand', 'pamid', 'varid', 'ref', 'alt', 'AF_info']]
        # to create tabix files, order must be chrom then pos
    
    # filter non-SNP
    _snp_letter = ('A', 'C', 'G', 'T', '.', -1)
    df = df.loc[(df['ref'].isin(_snp_letter)) & (df['alt'].isin(_snp_letter))]
    
    genome = Fasta(genome_path)
    ref_seqs, alt_seqs, ref_preds, alt_preds, abs_diff = [], [], [], [], []
    
    for i in range(len(df)): 
        chrom = df['chrom'].iloc[i] # exists for novar
        start = df['start'].iloc[i] # exists for novar
        end = df['end'].iloc[i] # exists for novar
        pos = df['pos'].iloc[i] # novar: == -1
        ref = df['ref'].iloc[i] # novar: == '.'
        alt = df['alt'].iloc[i] # novar: == '.'
        
        # run CROTON on ref_seq w/ or w/o variants
        seq = genome[chrom][start:end]
        if df['strand'].iloc[i] == '-':
            ref_seq = seq.reverse.complement.seq
        else: ref_seq = seq.seq
        ref_seqs.append(ref_seq)
        
        ref_oh = one_hot_encode(ref_seq, 'ACGT').reshape((1, 60, 4))
        ref_pred = model.predict(ref_oh).flatten()
        ref_preds.append(ref_pred)

        if (ref != '.'): # alternate variant cases (containing SNV)
            inx = pos - start
            assert seq.seq[inx-1] == ref
            
            # for initial analysis only focus on SNP
            # for more complex cases, use decan cases for handling
            if len(ref) > 1 or len(alt) > 1:
                continue
            
            # fill in alternative vars
            alt_seq = list(seq.seq)
            alt_seq[inx-1] = alt
            alt_seq = ''.join(alt_seq)

            if df['strand'].iloc[i] == '-': alt_seq = reverse_complement(alt_seq)
            alt_seqs.append(alt_seq)

            # make alternate prediction
            alt_oh = one_hot_encode(alt_seq, 'ACGT').reshape((1, 60, 4))
            alt_pred = model.predict(alt_oh).flatten()
            alt_preds.append(alt_pred)
            
            # get abs_diffs
            abs_diffs = (np.absolute(np.array(ref_pred) - np.array(alt_pred))) # list of absolute differences
            abs_diff.append(max(abs_diffs))
        
        else: # no variant (no SNV) cases
            alt_seqs.append('.')
            alt_preds.append(["{0:.0f}".format(-1)] * 6)
            abs_diff.append('.')
    
    # fill in df
    df['ref_seq'] = ref_seqs
    df['alt_seq'] = alt_seqs
    df['ref_preds'] = ref_preds
    df['alt_preds'] = alt_preds
    df['abs_diff'] = abs_diff
    
    tasks_ordering = ['del_frq', '1ins', '1del', 'onemod3', 'twomod3', 'frameshift']
    for i, t in enumerate(tasks_ordering):
        df['ref_%s'%t] = df['ref_preds'].apply(lambda x: x[i])
        df['alt_%s'%t] = df['alt_preds'].apply(lambda x: x[i])
    df = df.drop(['ref_preds', 'alt_preds'], axis=1)
    #df = df.sort_values(by='abs_diff', ascending=False)
    
    return df

def get_incomplete_bedfps():
    """
    Find files in bed_dir without a corresponding tsv file in tsv_dir
        - So that Slurm jobs do not run for files already created
    """
    bed_fp_lst = sorted([os.path.join(path, name) for path, subdirs, files in os.walk(bed_dir) for name in files])
    bed_gns = []
    for fp in bed_fp_lst:
        genename = fp.replace(bed_dir, '')
        genename = genename.replace('.bed', '')
        bed_gns.append(genename)
    
    tsv_fp_lst = sorted([os.path.join(path, name) for path, subdirs, files in os.walk(tsv_dir) for name in files]) # already created
    tsv_gns = []
    for fp in tsv_fp_lst:
        genename = fp.replace(tsv_dir, '')
        genename = genename.replace('.tsv', '')
        tsv_gns.append(genename)
    
    incomplete_gns = np.setdiff1d(bed_gns,tsv_gns)
    incomplete_bedfps = sorted([bed_dir + gn + '.bed' for gn in list(incomplete_gns)])
    print(len(incomplete_bedfps)) #787, 196, 24
    #print(incomplete_bedfps)
    return incomplete_bedfps


def get_tsvs(chrom, incomplete=False, model_fp=model_fp):
    """
    Run get_varscore function on every file in bed_dir
        - final .tsv dataframes have the columns
            ['chrom', 'pamid', 'genename', 'num', 'start', 'end', 'strand', 'pos',
            'varid', 'ref', 'alt', 'inx', 'ref_seq', 'alt_seq', 'ref_preds',
            'alt_preds', 'abs_diff', 'ref_del_frq', 'alt_del_frq', 'diff_del_frq',
            'ref_1ins', 'alt_1ins', 'diff_1ins', 'ref_1del', 'alt_1del',
            'diff_1del', 'ref_onemod3', 'alt_onemod3', 'diff_onemod3',
            'ref_twomod3', 'alt_twomod3', 'diff_twomod3', 'ref_frameshift',
            'alt_frameshift', 'diff_frameshift'] 
    
    Parameters 
    ----------
    chrom: ['X', 'Y', 1, 2 ... 22] 
        Chromosome to run on
    model_fp: str
        Filepath to model to make predictions with
    incomplete: bool
        If True--> will run get_incomplete_bedfps to get list of bed_fps 
            - Only tsvs that have not been made will be made
        If False--> will just run on all bed_fps detected in bed_dir
    
    Yield
    -----
    ([genename].tsv) files in (tsv_dir + /[chromosome]) directories
    
    Example
    -------
    >>> first couple of rows of PDCD1.tsv
      chrom pamid    genename  num      start        end strand        pos         varid ref alt  inx  ...  ref_1del  alt_1del diff_1del ref_onemod3  alt_onemod3  diff_onemod3  ref_twomod3  alt_twomod3  diff_twomod3  ref_frameshift  alt_frameshift  diff_frameshif0     chr2   PDCD1|57    PDCD1   57  241851253  241851313      +  241851283  rs1284638279   T   G   30  ...  0.032371  0.160016 -0.127645    0.699742     0.338150      0.361592     0.171198     0.383210     -0.212012        0.863511        0.718307         0.145204
      chr2  PDCD1|148   PDCD1  148  241852758  241852818      -  241852789   rs535799968   C   T   31  ...  0.135946  0.036195  0.099751    0.252422     0.464601     -0.212179     0.299638     0.151222      0.148416        0.556249        0.611154        -0.054905
      chr2  PDCD1|70    PDCD1   70  241852175  241852235      +  241852205   rs141119263   G   A   30  ...  0.090900  0.024949  0.065952    0.457806     0.748258     -0.290451     0.320004     0.145324      0.174680        0.768263        0.885476        -0.117213
    """
    chrom = str(chrom)
    model = load_model(model_fp)
    if incomplete: fp_lst = get_incomplete_bedfps() # NOT running for the first time
    else: fp_lst = sorted([os.path.join(path, name) for path, subdirs, files in os.walk(bed_dir+'/%s/'%chrom) for name in files]) #20143 NOT 18079, running for first time

    for fp in fp_lst:
        genename = fp.split('/')[-1].replace('.bed', '')
        print(genename)
        df = get_varscore(bed_fp=fp, model=model)
        df.to_csv(tsv_dir + '/%s/%s.tsv'%(chrom, genename), index=False, sep="\t")


if __name__ == "__main__": # see sbatch_tsv.sh
    parser = argparse.ArgumentParser(description='Making tsvs')
    parser.add_argument('--chrom', help='Chromosome with which to make bedsn')
    args = parser.parse_args()
    make_tsvdirs()
    get_tsvs(chrom=args.chrom)


# DID NOT FILTER DISRUPTIVE PAMs WHEN MAKING -- COULD BE HELPFUL DOWNSTREAM TO FIGURE OUT WHICH HAVE DISRUPTIVE PAMS
def filter_disruptive_pams(genename=None, bed_df=None): # get rid of SNVs that disrupt PAM site
    if bed_df is None:
        bed_df = pd.read_csv("./data/data/genome/%s.tsv" % genename, sep='\t')
    posdrop_series = (bed_df['strand']=='+') & ((bed_df['inx'] == 35)|(bed_df['inx'] == 36)) # interferes with PAM on + strand seqs
    bed_df = bed_df[~posdrop_series]
    negdrop_series = (bed_df['strand']=='-') & ((bed_df['inx'] == 25)|(bed_df['inx'] == 26)) # interferes with PAM on - strand seqs
    bed_df = bed_df[~negdrop_series]
    
    return bed_df

'''
AC253536.7
     chrom          pamid     start       end strand  pos varid ref alt AF_info
323  chr22  AC253536.7|12  23973735  23973795      +   -1     .   .   .      {}
324  chr22  AC253536.7|13  23973741  23973801      +   -1     .   .   .      {}
325  chr22  AC253536.7|14  23973747  23973807      -   -1     .   .   .      {}
326  chr22  AC253536.7|15  23973747  23973807      +   -1     .   .   .      {}
'''
