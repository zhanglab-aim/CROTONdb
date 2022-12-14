# Test get_CDSpams.py functions using pytest

# List of functions to test:
#  1. backend.get_CDS_df
#  2. backend.get_CDS_seq
#  3. backend.get_PAM_coords
#  4. backend.make_CDSpamsbed
#  5. backend.get_ref_PAM_seq
#  6. backend.subset_CDSpams

from get_CDSpams import getCDSpams

def test_get_CDS_df():
    #Test for: sequence end (end) - sequence start (start) = sequence length (strand)
    
    #Test for: Variable 'end' is greater than 'start'
    pass

def test_get_CDS_seq():
    #Test that only Chr1-22, chrX, CDS in pam_df

    #Test that PAM columns were added to pam_df.csv
    pass

def test_get_PAM_coords():
    #
    pass

def test_make_CDSpamsbed():
    pass

def test_get_ref_PAM_seq():
    pass

def test_subset_CDSpams():
    pass

