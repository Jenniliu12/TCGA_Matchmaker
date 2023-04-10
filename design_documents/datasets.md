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


3. Real dataset for answering a biological question using the tool
(size, structure, basic stats)

Gene expression profiles can be downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&type=2) or the TCGA website (https://portal.gdc.cancer.gov/repository). The datasets are likely to contain excess information. The correct data is to be filtered.
The data pre-processing and reformatting is to be conducted by the user prior to uploading to the software. The correct format is a comma-separated-file with the first column ("symbol") containing gene names and the second column ("value") containing the expression levels. 
For the firrst prototype, both profiles are to be uploaded by the user. In the future, the tool aims to have provided a database of representative cancer samples.
The average size of the profiles is only a few MBs containing 60,000 genes. The datasets are anticipated to be computationally light.
With two uploaded cancer samples, the question to be answered is how many genes are highly expressed in both cancers, indicating commonalities among different cancer types. 
