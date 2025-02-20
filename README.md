
# Reference management
Generate genome resources by downloading genome files and generating required files to facilitate running nf-core workflows.

## Description

This repository is to manage reference genomes. The scripts here are used to:

- Download reference genomes from public repositories (ENSEMBL and NCBI)
    - ENSEMBL: use API to get info for a given specie and download fasta and gtf files from the ftp server
    - NCBI: use datasets tool to get list of assemblies and download the fasta and gtf files
- Structure the reference files and link them to unique assembly names by write the paths to a config file that can be fed into nf-core/references workflow
- Generate star and bwa indices for each of the species using nf-core/references
- Write paths of the reference files and their indices to a main config file that can be used for data analysis with nf-core workflows.

### Get NCBI reference genomes:

Download all reference genomes from NCBI:

`python download_ncbi_data.py --search_term 'all' --download_species -m 1000000 --processes 40`
 
Search for a particular species and download it:

`python download_ncbi_data.py --search_term 'human' --download_species -m 1`
 

### Get ENSEMBL reference genomes

Get the list of all reference genomes in primates division from release 113 - no download:

`python download_ensembl_data.py --search_term all --division primates --release 113 --species_list_file ensembl_species_primates.txt`

Get the list of all reference genomes and Download all of them (fasta and gtf files) in primates division from release 113:

`python download_ensembl_data.py --download_species --search_term all --file_types fasta gtf --division primates --release 113 --processes 20`

-----

Get the list of all reference genomes in plants division from current release - no download:

`python download_ensembl_data.py --search_term all --division plants --release current --species_list_file ensembl_species_plants.txt`

Get the list of all reference genomes and Download all of them (fasta and gtf files) in plants division from the current release:

`python download_ensembl_data.py --download_species --search_term all --file_types fasta gtf --division plants --release current --processes 20`

-----

Download reference genomes (fasta and gtf files) for a given species (e.g. homo_sapiens) in primates division from release 113:

`python download_ensembl_data.py --search_term homo_sapiens --file_types fasta gtf --division primates --release 113`

Download reference genomes (fasta and gtf files) for a given species (e.g. actinidia_chinensis) in primates division from the release 60:

`python download_ensembl_data.py --search_term actinidia_chinensis --file_types fasta gtf --division plants --release 60`


Update the current release for all spcies in the plant division with the most recent files (force re-downloading all files): 

`python download_ensembl_data.py --search_term all --file_types fasta gtf --division plants --release current --download_species --force_replace`