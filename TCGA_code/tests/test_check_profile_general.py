from TCGA_code import match_computation as m
import pandas as pd
import numpy as np

# THIS TEST CASE IS A GENERAL CASE
def test_check_profile1():

    # Create two dataframes that represent the input reference profile and the TCGA data:
    gene_list = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5']
    value_list = [1, 4, 7, 12, 15]
    data_df = np.stack((gene_list, value_list), axis = 1)
    profile = pd.DataFrame(data_df, columns=["symbol", "value"])

    gene_list2 = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5']
    value_list2 = [1, 4, 7, 12, 15]
    data_df2 = np.stack((gene_list2, value_list2), axis = 1)   
    sample_data = pd.DataFrame(data_df2, columns=["symbol", "value"])

    # (1)
    # Test the function. We expect that all genes of the reference profile are present in the TCGA dataset.
    new_ref_profile, new_sample = m.check_profile(profile, sample_data, add_missing = False, output = False)
    # This assert statement is expected to be true as no genes should have been dropped here.
    assert new_ref_profile.iloc[:,0].tolist() == profile.iloc[:,0].tolist(), "All genes in the profile should be in the TCGA sample data!"

    # (2)
    # Test the function. We expect that all genes of the reference profile are present in the TCGA dataset.
    new_ref_profile, new_sample = m.check_profile(profile, sample_data, add_missing = True, output = False)
    # This assert statement is expected to be true as no genes should have been dropped here.
    assert new_sample.iloc[:,0].tolist() == sample_data.iloc[:,0].tolist(), "All genes in the profile should be in the TCGA sample data!"
    
def test_check_profile2():

    # Create two dataframes that represent the input reference profile and the TCGA data:
    gene_list = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5']
    value_list = [1, 4, 7, 12, 15]
    data_df = np.stack((gene_list, value_list), axis = 1)
    profile = pd.DataFrame(data_df, columns=["symbol", "value"])

    gene_list2 = ['gene1']
    value_list2 = [1]
    data_df2 = np.stack((gene_list2, value_list2), axis = 1)   
    sample_data = pd.DataFrame(data_df2, columns=["symbol", "value"])

    # (1)
    # Test the function. We expect that genes of the reference profile are missing in the TCGA dataset.
    # Here, the genes in the reference profile are expected to be dropped. That means we end with only 'gene1'
    new_ref_profile, new_sample = m.check_profile(profile, sample_data, add_missing = False, output = False)
    print(new_ref_profile)
    # This assert statement is expected to be true as genes should have been dropped here and both now contain only a single gene.
    #assert new_ref_profile.iloc[:,0].tolist() == profile.iloc[:,0].tolist(), "All genes in the profile should be in the TCGA sample data!"

    # (2)
    # Test the function. We expect that genes of the reference profile are missing in the TCGA dataset.
    new_ref_profile, new_sample = m.check_profile(profile, sample_data, add_missing = True, output = False)
    # This assert statement is expected to be true as genes have been added here.
    #assert new_sample.iloc[:,0].tolist() == sample_data.iloc[:,0].tolist(), "All genes in the profile should be in the TCGA sample data!"

    



