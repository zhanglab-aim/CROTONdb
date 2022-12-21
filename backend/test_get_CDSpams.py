# Test get_CDSpams.py functions using pytest

# List of functions to test:
#  1. backend.get_CDS_df
#  2. backend.get_CDS_seq
#  3. backend.get_PAM_coords
#  4. backend.make_CDSpamsbed
#  5. backend.get_ref_PAM_seq
#  6. backend.subset_CDSpams

import pytest
from . import configs
from .get_CDSpams import get_CDS_df, get_CDS_seq, get_PAM_coords

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
    #
    pass

#def test_make_CDSpamsbed():
#    pass

#def test_get_ref_PAM_seq():
#    pass

#def test_subset_CDSpams():
#    pass

