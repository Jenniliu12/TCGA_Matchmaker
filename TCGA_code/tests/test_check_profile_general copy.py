from TCGA_code import match_computation as m
import pandas as pd

# THIS TEST CASE IS A GENERAL CASE
def test_check_profile_gen():
    profile = ['ABCD', 'DEFG', 'QWERT', 'QWERT', 'ASDFG', 'PPPPP']
    sample_data = ['ABCD', 'DEFG', 'QWERT', 'ASDFG', 'PPPPP']

    is_present = m.check_profile(profile, sample_data)
    truth = False
    assert is_present == truth, "All genes in the profile should be in the sample data!"

# THIS TEST CASE IS A CASE, WHEN NOT ALL THE SAMPLE_DATA IS IN THE PROFILE
def test_check_profile_gen2():
    profile = ['ABCD', 'DEFG', 'QWERT', 'QWERT', 'ASDFG', 'OOOII' ]
    sample_data = ['ABCD', 'DEFG', 'QWERT', 'ASDFG', 'PPPPP']

    is_present = m.check_profile(profile, sample_data)
    truth = False
    assert is_present == truth, "Not all genes in the profile should be in the sample data!"


# THIS TEST CASE IS A CASE, WHEN THE SAMPLE DATA IS EMPTY, which should output FALSE, and print a statement. 
def test_check_profile_gen3():
    profile = ['ABCD', 'DEFG', 'QWERT', 'QWERT', 'ASDFG', 'OOOII' ]
    sample_data = []

    is_present = m.check_profile(profile, sample_data)
    truth = True
    assert is_present == truth, "All genes in the profile should be in the sample data!"