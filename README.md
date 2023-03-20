# TCGA_Matchmaker_2
Version 2 of the BIOINF575 Project: A tool to match gene expression profiles to find cancerous biomarkers.

reference to example data that can be used to run your package and 
what the expected result would be and 
write instruction on how a user will run the package on that data add this to the README.md file

An example pipeline of this algorithm is provided in the 'TCGA_code/example_test' folder:
The necessary and already formatted input files are provided. 
"example_input_sample.txt" represents the query expression profile and "example_profile.txt" represents the comprehensive TCGA cancer profiles.
Their formats need to be in the following format: two tab-delimeted columns. Gene_id/name on the left; expression levels on the right.

The example pipeline simulates the main components of the software:
Reading data in. 
The expected result of the distance computation is 0.319. 

