# Description of steps required for the Snakemake demo for mapping
# Please read this file to the end before you start running commands so you know what to expect. 
# Marked with >>> are commands that you will run in the command line

# The data folder contains the reference genome and index files needed for alignment and a set of reads in the sample folder. 

# Open the command line and navigate to the location of this README.txt file.

# Create the environment for the workflow:
# The file envs/mapping.yaml contains the configuration for what the environment needed for this workflow to run will contain and where to get the required tools from. The file environment.yaml contains a settings for an environment with more tools.

# mapping.yaml 
# channels:
#  - bioconda
#  - conda-forge
#  - defaults
# dependencies:
#  - snakemake-minimal =5.2.4
#  - bwa =0.7.17
#  - samtools =1.9


# The command that follows will create a folder snakemake-mapping on your anaconda installation directory, for me the folder was at: /Users/mitrea/anaconda3/envs/snakemake-mapping.
# If you want to start over and rerun the command above look for that folder in your system and remove it.


>>> conda env create --name snakemake-mapping --file envs/mapping.yaml


# In order to use the tools snakemake, bwa and samtools  from the snakemake-mapping environment you need to activate it by running the following command.


>>> conda activate snakemake-mapping


# To deactivate the environment you can use: conda deactivate

# Navigate to the align_reads_example1 folder


>>> cd align_reads_example1


# To run the alignment and produce a .bam file you tell the workflow what output to create and it will work its way back to the input, run the rules that are set in sequence and create the output.
# It will run the commands in the shell part of the rule here alignment of single-end reads with bwa and then the output is converted to a .bam file (-b option) with samtools view and saved in the mapped_reads folder.
# The rule(s) in the Snakemake file describe the workflow.  
# It will run the commands in the shell part of the rules on the files in the input part of the rule and create the files in the output part of the rules. 

# Snakefile

# rule bwa_map:
#    input:
#        "../data/genome.fa",
#        "../data/samples/{sample}.fastq"
#    output:
#        "../mapped_reads/{sample}.bam"
#    shell
#        "bwa mem {input} | samtools view -b - > {output}"


>>> snakemake ../mapped_reads/A.bam


# You can provide multiple files and the workflow will apply the rules to each input


>>> snakemake ../mapped_reads/{B,C}.bam


# If no error, and you see the files A.bam, B.bam, C.bam in the mapped_reads folder you have successfully ran the workflow!
# Congratulations! 

# You can create new workflows that follow the same structure.
# I recommend creating the workflow (Snakefile) in a new folder because the file that the tool snakemake uses has to be called Snakefile.

# An example of a workflow with 2 steps is provided in the align_reads_example2 folder.
# The Snakefile provides the rules of the workflow.
# The rule map_reads - does alignment of paired-end reads with bwa and then the output is converted to a .bam file (-b option) with samtools view and saved in the mapped folder this time
# The rule bam_view - converts the input file to a .sam file with samtools view. Sam files are plain text files that you can view with a text editor.
# https://samtools.github.io/hts-specs/SAMv1.pdf

# Snakefile

# rule map_reads:
#     input:
#         "../data/genome.fa",
#         "../data/samples/{sample}_1.fq",
#         "../data/samples/{sample}_2.fq"
#     output:
#         "../mapped/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -b - > {output}"
# 
# rule bam_view:
# 	input:
# 		"../mapped/{sample}.bam"
# 	output:
# 		"../mapped/{sample}.txt"
# 	shell:
# 		"samtools view {input} > {output}"


# Navigate to the align_reads_example2 folder and run the new workflow on the paired-end reads for sample D.


>>> cd ../align_reads_example2


>>> snakemake ../mapped/D.txt


>>> cd ..


# If no error, and you see the D.bam and D.txt files in the mapped folder you have successfully ran the workflow!
# Congratulations! 

DONE!


The commands from this file are also available in the commands_to_run.sh which you can run as a bash script.


https://snakemake.readthedocs.io/en/stable/
https://snakemake.readthedocs.io/en/v3.10.0/tutorial/welcome.html
https://snakemake.readthedocs.io/en/v3.10.1/tutorial/basics.html
https://snakemake.bitbucket.io/snakemake-tutorial.html




