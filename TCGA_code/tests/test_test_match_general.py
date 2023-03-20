from TCGA_code import match_computation as m
import pandas as pd

def test_test_match_gen():
    profile = ['ABCD', 'DEFG', 'QWERT', 'QWERT', 'ASDFG', 'PPPPP']
    sample_data = ['ABCD', 'DEFG', 'QWERT', 'ASDFG', 'PPPPP']

    threshold = 1

    is_match = m.test_match(profile, sample_data, threshold)
    truth = False

    assert is_match == truth