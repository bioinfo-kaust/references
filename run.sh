#!/bin/bash
#SBATCH -p batch
#SBATCH --time=90:00:00
#SBATCH --mem=100GB
#SBATCH --ntasks=20

# 1. Download Ensembl genome files

# 2. Download NCBI genome files
python download_ncbi_data.py --download_species --accepted_divisions Primates -m 100000 --file_types genome gtf --processes 40

# 3. Generate YAML file based on the genomes directory to list existing files per species

# 4. Run nf-core/references workflow to generate indices, intervals, etc

#5 . Generate genomes.config file to list the files for each species for running other nf-core workflows


