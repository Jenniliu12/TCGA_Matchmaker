from TCGA_code import match_computation as m
import pandas as pd

def test_compute_dist_gen():
    # profile = pd.Series([2,3,4], index = ["g1", "g2", "g3"])
    # sample_data = pd.Series([2,4,4,5], index = ["g1", "g2", "g4", "g5"])

    profile = pd.Series([2,3,4,7,8,9,1,4,6,1], index = ["g1", "g2", "g3", "g4", "g5", "g6","g7", "g8", "g9", "g10"])
    sample_data = pd.Series([2,4,4,58,9,4,2,8,6,1], index = ["g1", "g2", "g4", "g5", "g12", "g21", "g47", "g8", "g9", "g10"])

    expected_distance = 0.646 #0.6463447741067103
    distance = m.compute_distance(profile, sample_data)
    # print(distance, expected_distance)
    assert distance == expected_distance, "The distance between is 0.646 not " + str(distance)# "The distance between " + str(profile) + " and " + str(sample_data) +" is 0 not " + str(distance)