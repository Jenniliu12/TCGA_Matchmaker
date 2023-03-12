# Brief description of the datasets I will use for this project

1. TCGA data
We will retrieve TCGA gene expression data from: 
(link: GCD: Isoform quantification ?)
https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.disease_type%22%2C%22value%22%3A%5B%22ductal%20and%20lobular%20neoplasms%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22breast%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22genes.is_cancer_gene_census%22%2C%22value%22%3A%5B%22true%22%5D%7D%2C%22op%22%3A%22in%22%7D%5D%7D
Data might be for isoforms not genes so may need to summarize


2. A data file with gene expression profile
It will contain a matrix with gene symbol or id and gene expression on each row
Ideally it will be Gene NCBI entry ID for the gene ID and TPM values for the gene expression
This file will be provided by the user, I will have data for testing
