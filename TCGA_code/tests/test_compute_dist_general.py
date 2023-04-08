from TCGA_code import match_computation as m
import pandas as pd

def test_compute_dist_gen():

    profile = pd.DataFrame([['gene1', 1], ['gene2', 5], ['gene3', 8]], columns=["symbol", "value"])
    
    sample_data = pd.DataFrame([['gene1', 1], ['gene2', 5], ['gene3', 8]], columns=["symbol", "value"])

    expected_distance = 1.0

    distance = m.compute_distance(profile, sample_data)

    assert distance == expected_distance



