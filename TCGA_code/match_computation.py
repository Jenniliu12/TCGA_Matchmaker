import numpy as np
import pandas as pd
import bamnostic as bs


def read_expr_profile(file_name):
	'''
	Read in the gene expression profile to be analyzed.

	Parameters: 
		file_name (string): the file name of the gene expression profile input file. (Most likely a gct file instead of a csv ?)
							Has to be a two column file. 

	Returns:
		Reference profile data (panda series): Series where labels are gene symbols and values are the expression of the respective gene

		
	d = {'a': 1, 'b': 2, 'c': 3}
	ser = pd.Series(data=d, index=['a', 'b', 'c'])
	-->
	a    1
	b    2
	c    3

	'''

	# Read input file:
	profile_df = pd.read_csv(file_name, sep = "\t", squeeze = True)

	# index = gene symbols
	index = ref_profile_data.iloc[:,1]

	# values = gene_expression
	values = ref_profile_data.iloc[:,2]

	ref_profile_data = pd.series(index, values)
	return ref_profile_data


def read_TCGA_sample(file_name):
	'''
	Read in the TCGA sample to be analyzed.

	Parameters: 
		file_name (string): the file name of the input file. A .bam file.

	Returns:
		sample_data (list of strings): a string of the query RNA sequence.
	
	'''

	# Read input file:
	profile_df = pd.read_csv(file_name, sep = "\t", squeeze = True)

	# index = gene symbols
	index = profile_df.iloc[:,1]

	# values = gene_expression
	values = profile_df.iloc[:,2]

	sample_profile_data = pd.series(index, values)
	return sample_profile_data 
	

def test_match(profile, sample_data, threshold):

	'''
	Checks the match of the sample data and the profile given a certain threshold.

	Parameters: 
		profile (list of strings): A list of genes from the profile
		sample_data (list of strings): A list of genes from the sample_data
		threshold (float): Number that determines whether or not the match is strong enough to be an actual match.

	Returns:
		is_match (bool): Wether or not all genes from the sample data are also present in the profile
	
	'''
	is_match = False
	is_present = check_profile(profile, sample_data)
	if is_present:
		distance = compute_distance(profile, sample_data)
		is_match = distance < threshold
	return is_match


def check_profile(profile, sample_data):
	# is_present = False
	# check if all genes in the profile are present

	'''
	Check if all genes in the profile are present

	Parameters: 
		profile (list of strings): A list of genes from the profile
		sample_data (list of strings): A list of genes from the sample_data

	Returns:
		is_present (bool): Whether or not all genes from the sample data are also present in the profile
	'''
	# Turn the lists into sets to remove duplicates. IS THIS NECESSARY SINCE THE GENE SYMBOLS ARE THE INDICES IN THE pd.series???
	profile = list(set(profile))
	sample_data= list(set(sample_data))

	# Check if one of the datasets is empty. s
	if len(sample_data) or len(profile) == 0:
		print('At least one of the data sets is empty.')
		is_present = False
		return is_present
	else:
		for data in sample_data: 
			if data not in profile:
				is_present = False
				return is_present

		is_present = True
		return is_present


def compute_distance(profile, sample_data):

	'''
	Compute the distance between the reference expression profile and a sample expression profile

	Parameters: 
		profile (pandas Series): a pandas series with the reference profil
								data is the expression level, rownames are gene symbols or IDs
		sample_data(pandas.Series): a pandas series with the sample profile data is the expression level, rownames are gene symbols or IDs

	Returns (float): match score - a distance type metric that shows how 
					0 would be minimum and 1 would be maximum
	
	'''
	distance = 1
	# compute distance here
	#distance = sample_data - profile
	set1 = set(sample_data.index)  # results in removal of duplicates
	set2 = set(profile.index)      # results in removal of duplicates
	intersection = list(set1.intersection(set2))  # returns the elements that are the same among set1 and set2
	overlap_no = len(intersection)  # returns the size of the intersection

	smaller_set_no = min(len(set1), len(set2))  # len = 3
	perc_overlap = overlap_no/smaller_set_no  # 1/3 = 0.33

	#distance = sample_data[intersection].corr(profile[intersection])
	#distance = round(distance, 3)
	
	if (perc_overlap > 0.6):
		print("There is a significant overlap!")
		#distance = sum(abs(sample_data[intersection]-profile[intersection]))  # TAKE CORRELATION
		distance = sample_data[intersection].corr(profile[intersection])
		distance = round(distance, 3)
	return distance
	