import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr   
from scipy.stats import zscore
import matplotlib.pyplot as plt


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
    assert len(ref_profile.axes[1]) == 2, "The number of columns must be 2: 'symbol, value'."
    assert np.isnan(any(ref_profile.iloc[:,0])) == False, "There are not defined gene names in the reference profile."
    assert np.isnan(any(ref_profile.iloc[:,1])) == False, "There are not defined expression levels in the reference profile."

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
    assert len(TCGA_profile.axes[1]) == 2, "The number of columns must be 2: 'symbol, value'"
    assert np.isnan(any(TCGA_profile.iloc[:,0])) == False, "There are not defined gene names in the TCGA profile."
    assert np.isnan(any(TCGA_profile.iloc[:,1])) == False, "There are not defined expression levels in the TCGA profile."

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
        profile.iloc[:,1] =(profile.iloc[:,1]-profile.iloc[:,1].min())/(profile.iloc[:,1].max()-profile.iloc[:,1].min())
        return profile
    elif method == "raw":
        return profile
    else:
        print("Not a valid input method")
        return 0
    
def expression_analysis(profile, sample_data, sensitivity_threshold = 0.05):
    '''
    Analyses similiarities between the reference expression profile and a sample expression profile:
        -> Top similiar genes within a threshold
    
    Input profiles are to be used after check_profile and check_TCGA
    Input profile gene columns are expected to match.

    Parameters: 
        profile (pd.DataFrame): Reference Profile dataset
        sample_data (pd.DataFrame): TCGA dataset

    Returns (float): 
        gene_ratio (pd.DataFrame) 
        similiar (list of genes that have similiar expression levels)

    ''' 

    assert len(profile.axes[0]) == len(sample_data.axes[0]), "Input datasets are not of same length."
    assert all(profile.eq(sample_data, 0)) == True, "Input gene profiles do not match. Check if they are sorted. Check if genes are missing."
    
    # Loop through the dfs and find the genes whose expression levels are most similiar within a certain threshold
    similiar = []  # list of genes that have similiar expression levels among two profiles
    ratio_list = []  # list ofthe ratio of genes  
    for gene in range(len(profile.axes[0])):
        ratio = profile.iloc[gene, 1]/sample_data.iloc[gene, 1]
        # Consider the nans. You have to edit this later! I want it to be representative. 
        if np.isnan(ratio):
            ratio = 0
        ratio_list.append(ratio)
        if abs(ratio-1) <= sensitivity_threshold:
            similiar.append(gene)

    gene_ratio = pd.DataFrame(list(zip(profile.iloc[:,0].tolist(), ratio_list)), columns =['symbol', 'ratio'])

    return gene_ratio, similiar


def gene_bar_chart(gene_ratio, show = "percentage"):
    '''
    Prints a bar chart with gene ratios. 
    
    Parameters: 
        gene_ratio (pd.DataFrame) = obtained from expression_analysis function.
        show (string) = determines whether percentage or discrete count values are shown

    Returns (float): 
        na
        (prints bar chart)

    ''' 
    # loop through dataframe to create bins
    category = []
    for ratio in gene_ratio["ratio"]:
        if ratio <= 0.2:
            category.append("0.0-0.2")
        elif ratio > 0.2 and ratio <= 0.4:
            category.append("0.2-0.4")
        elif ratio > 0.4 and ratio <= 0.6:
            category.append("0.4-0.6")
        elif ratio > 0.6 and ratio <= 0.8:
            category.append("0.6-0.8")
        elif ratio > 0.8 and ratio <= 1:
            category.append("0.8-1.0")
        elif ratio > 1.0 and ratio <= 1.2:
            category.append("1.0-1.2")
        elif ratio > 1.2 and ratio <= 1.4:
            category.append("1.2-1.4")
        elif ratio > 1.4 and ratio <= 1.6:
            category.append("1.4-1.6")
        elif ratio > 1.6 and ratio <= 1.8:
            category.append("1.6-1.8")
        elif ratio > 1.8 and ratio <= 2.0:
            category.append("1.8-2.0")
        elif ratio > 2.0:
            category.append(">2.0")
        else:
            category.append("nan")
    
    # Count the values of each bin:
    gene_ratio["category"] = category
    gene_ratio_count = gene_ratio.groupby('category')["category"].count()


    p1 = plt.bar(gene_ratio_count.index.sort_values(), gene_ratio_count,  color=['#FFA500', '#FFA500', '#FFA500', '#FFA500', 
                                                                                 '#2E8B57', '#2E8B57', 
                                                                                 '#191970', '#191970','#191970', '#191970','#191970','#191970'])
                                                                                 

    for rect1 in p1:
        height = rect1.get_height()
        if show == "count":
            plt.annotate("{}".format(height),(rect1.get_x() + rect1.get_width()/2, height+.05),ha="center",va="bottom",fontsize=9)
        else: # shows percentage
            plt.annotate("{}%".format(round(height/len(gene_ratio),3)),(rect1.get_x() + rect1.get_width()/2, height+.05),ha="center",va="bottom",fontsize=9)

    plt.xticks(rotation="vertical")
    plt.show()






            
