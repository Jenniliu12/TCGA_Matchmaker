#!/bin/bash

conda env create --name snakemake-mapping --file envs/mapping.yaml

conda activate snakemake-mapping

cd align_reads_example1

snakemake ../mapped_reads/A.bam

snakemake ../mapped_reads/{B,C}.bam

cd ../align_reads_example2

snakemake ../mapped/D.txt

cd ..
