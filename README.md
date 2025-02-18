# Reference management

Generate genome resources by downloading genome files and generating required files to facilitate running nf-core workflows.

## Description

This repository is to manage reference genomes. The scripts here are used to:
- Download reference genomes from public repositories (ENSEMBL and NCBI)
  - ENSEMBL: use API to get info for a given specie and download fasta and gtf files from the ftp server
  - NCBI: use `datasets` tool to get list of assemblies and download the fasta and gtf files
- Structure the reference files and link them to unique assembly names by write the paths to a config file that can be fed into nf-core/references workflow
- Generate star and bwa indices for each of the species using nf-core/references
- Write paths of the reference files and their indices to a main config file that can be used for data analysis with nf-core workflows.

## Further information

Visit our wiki page for further information : https://bclwiki.kaust.edu.sa/en/bix/analysis/public/references


**Contact:**
- Husen Umer for reference genomes

