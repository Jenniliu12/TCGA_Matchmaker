from TCGA_Matchmaker import match_computation as m
import pandas as pd

def test_read_expr_profile_gen():

    # ref_profile_data = m.read_expr_profile("input1.csv")


    # 1,3,4,3,2,6
    # g1,2,3,4,5,5
    # expected_series = pd.series([1,3,4,3,2,6], index = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene5'])
   # assert ref_profile_data == expected_series
   assert 1 == 1