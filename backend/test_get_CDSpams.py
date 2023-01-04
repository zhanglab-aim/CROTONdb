# Test get_CDSpams.py functions using pytest

# List of functions to test:
#  1. backend.get_CDS_df
#  2. backend.get_CDS_seq
#  3. backend.get_PAM_coords
#  4. backend.make_CDSpamsbed
#  5. backend.get_ref_PAM_seq
#  6. backend.subset_CDSpams

import pandas as pd
import pytest
from . import configs
from .get_CDSpams import get_CDS_df, get_CDS_seq, get_PAM_coords, make_CDSpamsbed

def test_get_CDS_df():
    cds_df = get_CDS_df(db_path=configs.GTF_DB_PATH)
    #Test that sequence end (end) - sequence start (start) is a positive integer
    CDS_df_rows=cds_df.shape[0]
    for x in range(CDS_df_rows):
        print ("row: ",x)
        assert ((((cds_df.iloc[x]['end'])-(cds_df.iloc[x]['start'])))>0)

    return cds_df


def test_get_CDS_seq(cds_df=test_get_CDS_df()):
    #Test that the sequence added to df is the correct length based on start and end locations using 
    cds_df = get_CDS_seq(df=cds_df, genome_path=configs.GENOME_FA_PATH)

    CDS_df_rows=cds_df.shape[0]
    for x in range(CDS_df_rows):
        assert ((((cds_df.iloc[x]['end'])-(cds_df.iloc[x]['start']))) == len(cds_df.iloc[x]['seq'])) 

    return cds_df

def test_get_PAM_coords(cds_df=test_get_CDS_seq()):

    #Test if PAM coordinates point to NGG nucleotides
    pam_df = get_PAM_coords(df=cds_df)

    for x in range (len(pam_df['pams'].iloc[0])):
        pam_loc=pam_df['pams'].iloc[0][x]
        pam=pam_df['seq'].iloc[0][pam_loc:pam_loc+3]
        assert pam == 'AGG' or 'CGG' or 'TGG' or 'GGG'

    return pam_df


def test_make_CDSpamsbed():

    #Test if all sequences have 60bp
    #Note: Should potential add an error code to detect when the start location is negative or when
    #the end location is > the sequence length. In either scenario, it indicates that the PAM is too close
    #either the start or end of the sequence in order 
    CDSpamsbed_df, skipped_genes = make_CDSpamsbed(pam_df=test_get_PAM_coords())

    print(CDSpamsbed_df)

    for x in range (len(CDSpamsbed_df)):
        seq_len=(CDSpamsbed_df['end'].iloc[x]) - (CDSpamsbed_df['start'].iloc[x])
        assert seq_len == 60

#def test_get_ref_PAM_seq():

    #Test if start and end sequence locations still match the reference sequence in length and content
#    pass

#def test_subset_CDSpams():
    #Test if bed file contents match the appropriate chromosomes
#    pass

