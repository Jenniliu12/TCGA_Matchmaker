import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr   
from scipy.stats import zscore


def read_expr_profile(file_name):
    '''
    Read in the gene expression profile to be analyzed. This GE profile is analyzed against a dataset of cancer samples.

    Removes duplicates and average these duplicate values.
    Asserts that the file is in the correct format of 2 columns. 
    Asserts that there are no undefined values in the value column

    Parameters: 
        file_name (string): the file name of the gene expression profile input file. (.csv file or path to .csv file)
                             Has to be a two column file in the format: 'symbol' (string), 'value' (float or integer)

    Returns:
        Reference profile data (panda df): Pandas dataframe with gene symbol and gene expression levels (pre-processed)

    '''
    
    ref_profile = pd.read_csv(file_name, sep = ";", encoding="UTF-8")

    # Assert that the column size is 2.
    assert ref_profile.axes[1] == 2, "The number of columns must be 2: 'symbol, value'."
    assert np.isnan(ref_profile.iloc[:,0]) == False, "There are not defined gene names in the reference profile."
    assert np.isnan(ref_profile.iloc[:,1]) == False, "There are not defined expression levels in the reference profile."

    # Average the duplicate values
    ref_profile = ref_profile.groupby('symbol', as_index=False).mean()

    return ref_profile


def read_TCGA_sample(file_name):
    '''
    Read in the TCGA profile to be analyzed. This TCGA profile serves as a dataset of cancer samples.

    Removes duplicates and average these duplicate values.
    Asserts that the file is in the correct format of 2 columns. 
    Asserts that there are no undefined values in the value column

    Parameters: 
        file_name (string): the file name of the TCGA input file. (.csv file or path to .csv file)
                             Has to be a two column file in the format: 'symbol' (string), 'value' (float or integer)

    Returns:
        Reference profile data (panda df): Pandas dataframe with gene symbol and gene expression levels (pre-processed)

    '''

    TCGA_profile = pd.read_csv(file_name, sep = ";")

    # Assert that the column size is 2.
    assert TCGA_profile.axes[1] == 2, "The number of columns must be 2: 'symbol, value'"
    assert np.isnan(TCGA_profile.iloc[:,0]) == False, "There are not defined gene names in the TCGA profile."
    assert np.isnan(TCGA_profile.iloc[:,1]) == False, "There are not defined expression levels in the TCGA profile."

    # Average the duplicate values
    TCGA_profile = TCGA_profile.groupby('symbol', as_index=False).mean()

    return TCGA_profile 




# def test_match(profile, sample_data, threshold):

#     '''
#     Checks the match of the sample data and the profile given a certain threshold.

#     Parameters: 
#         profile (list of strings): A list of genes from the profile
#         sample_data (list of strings): A list of genes from the sample_data
#         threshold (float): Number that determines whether or not the match is strong enough to be an actual match.

#     Returns:
#         is_match (bool): Wether or not all genes from the sample data are also present in the profile

#         I don´t actually think I need this?

#     '''
#     is_match = False
#     is_present = check_profile(profile, sample_data)
#     if is_present:
#         distance = compute_distance(profile, sample_data)
#         is_match = distance < threshold
#     return is_match


def check_profile(profile, sample_data, add_missing = False, output = True):
    '''
    Check if all genes in the reference profile are present in the TCGA dataset. 
    If there are genes missing in the TCGA dataset (that are present in the ref-profile data):
        A: Drop the genes from the reference profile data.
        B: Add the genes to the TCGA dataset with zero expression values.

    Parameters: 
        profile (pd.DataFrame): Reference Profile dataset
        sample_data (pd.DataFrame): TCGA dataset
        add_missing (boolean): Signifies drop or add of data.
        output (boolean): Outputs the changes happening to the dataset.

    Returns:
    
        profile (pd.DataFrame): Updated profile. 
                                If add_missing = True: Unchanged
                                If add_missing = False: Genes are dropped from ref-profile.
                                
        sample_data (pd.DataFrame): Updated profile.
                                If add_missing = True: New genes with zero expression levels are added.
                                If add_missing = False: Unchanged.
    '''

    # CHECK IF ALL GENES IN THE PROFILE ARE PRESENT IN THE TCGA SAMPLE DATA.
    
    # profile and sample_data are dfs. Extract their gene names as lists.
    profile_genes = profile.iloc[:,0].tolist()
    TCGA_genes = sample_data.iloc[:,0].tolist()

    # Check if all values in profile_genes are also in TCGA_genes
    missing_genes = []
    missing_genes_idx = []
    idx = 0
    for gene in profile_genes:
        if gene not in TCGA_genes:
            missing_genes.append(gene)
            missing_genes_idx.append(idx)
        idx += 1 
    
    # If we don´t want to add temporary zero expression levels to the TCGA_dataset, just remove them from
    # the reference profile set and report them as missing.
    add_genes_dict = {}
    if add_missing == False:
        # Drop the missing genes from the reference profile
        profile = profile.drop(missing_genes_idx)
        
        if output == True:
            print("The following genes are dropped from the reference profile dataset and not considered as they do not exist in this TCGA sample database:")
            print(missing_genes)
        
        return profile, sample_data, missing_genes
    else:
        for gene in missing_genes:
            add_genes_dict[gene] = 0
            

        add_genes = pd.DataFrame(add_genes_dict.items(), columns=['symbol', 'value'])
        
        # Append the missing genes with 0 expression levels and sort genes alphabetically 
        sample_data = pd.concat([sample_data, add_genes], ignore_index=True) 
        sample_data = sample_data.groupby('symbol', as_index=False).mean()
        
        if output == True:
            print("Gene(s) from the reference profile are missing in the TCGA dataset and are added to the TCGA dataset with zero expression levels.")
            print(missing_genes)
        
        return profile, sample_data, missing_genes
        
        
def check_TCGA(profile, sample_data, add_missing = False, output = True):
    '''
    Check if all genes in the TCGA dataset are present in the reference profile.
    If there are genes missing in the reference profile dataset (that are present in the TCGA dataset):
        A: Drop the genes from the TCGA data.
        B: Add the genes to the reference profile dataset with zero expression values.

    Parameters: 
        profile (pd.DataFrame): Reference Profile dataset
        sample_data (pd.DataFrame): TCGA dataset
        add_missing (boolean): Signifies drop or add of data.
        output (boolean): Outputs the changes happening to the dataset.

    Returns:
    
        profile (pd.DataFrame): Updated profile. 
                                If add_missing = True: New genes with zero expression levels are added.
                                If add_missing = False: Unchanged.
                                
        sample_data (pd.DataFrame): Updated profile.
                                If add_missing = True: Unchanged.
                                If add_missing = False: Genes are dropped from TCGA dataset.
    '''
    
    # CHECK IF ALL GENES IN THE TCGA ARE PRESENT IN THE REFERENCE PROFILE DATASET.
    
    # profile and sample_data are dfs. Extract their gene names as lists.
    profile_genes = profile.iloc[:,0].tolist()
    TCGA_genes = sample_data.iloc[:,0].tolist()

    # Check if all values in TCGA_genes are alson in profile_genes
    missing_genes = []
    missing_genes_idx = []
    idx = 0
    for gene in TCGA_genes:
        if gene not in profile_genes:
            missing_genes.append(gene)
            missing_genes_idx.append(idx)
        idx += 1

        
    # If we don´t want to add temporary zero expression levels to the ref_profile_dataset, just remove them from
    # the TCGA profile set and report them as missing.
    add_genes_dict = {}
    if add_missing == False:
        # Drop the missing genes from the reference profile
        sample_data = sample_data.drop(missing_genes_idx)
        
        if output == True:
            print("The following genes are dropped from the TCGA dataset and not considered as they do not exist in this Reference profile dataset:")
            print(missing_genes)
        
        return profile, sample_data, missing_genes
    else:
        for gene in missing_genes:
            add_genes_dict[gene] = 0
            

        add_genes = pd.DataFrame(add_genes_dict.items(), columns=['symbol', 'value'])
        
        # Append the missing genes with 0 expression levels and sort genes alphabetically 
        profile = pd.concat([profile, add_genes], ignore_index=True) 
        profile = profile.groupby('symbol', as_index=False).mean()
        
        if output == True:
            print("Gene(s) from the TCGA profile are missing in the reference dataset and are added to the reference dataset with zero expression levels.")
            print(missing_genes)

        return profile, sample_data, missing_genes




def compute_distance(profile, sample_data):

    '''
    Compute the distance between the reference expression profile and a sample expression profile as a pearson correlation

    Parameters: 
        profile (pandas df): a pandas df with the reference profil
                                data is the expression level, rownames are gene symbols or IDs
        sample_data(pandas df): a pandas df with the sample profile data is the expression level, rownames are gene symbols or IDs

    Returns (float): match score - correlation of expression values distance type metric that shows how similiar the gene expression data sets are:
                    0 would be minimum and 1 would be maximum

    '''
    profile_levels = profile.iloc[:,1].tolist()
    TCGA_levels = sample_data.iloc[:,1].tolist()
    
    assert len(profile_levels) == len(TCGA_levels), "Input datasets are not of same length."
    
    distance, p = pearsonr(profile_levels, TCGA_levels)
    return round(distance,4)


def normalize_profile(profile, method):

    # SOURCE for mean and min-max: https://stackoverflow.com/questions/26414913/normalize-columns-of-a-dataframe
    if method == "z-score":
        profile.iloc[:,1] = zscore(profile.iloc[:,1])
        return profile
    
    elif method == "mean":
        profile.iloc[:,1] = (profile.iloc[:,1]-profile.iloc[:,1].mean())/profile.iloc[:,1].std()
        return profile
    
    elif method == "min-max":
        profile.iloc[:,1] =( profile.iloc[:,1]-profile.iloc[:,1].min())/(profile.iloc[:,1].max()-profile.iloc[:,1].min())
        return profile
    
    else:
        print("Not a valid input method")
        return profile