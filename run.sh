#!/bin/bash
#SBATCH -p batch
#SBATCH --time=90:00:00
#SBATCH --mem=200GB
#SBATCH --ntasks=40

module load nextflow
module load singularity
#make sure datasets binary from NCBI is installed

# 1. Download Ensembl genome files
python download_ensembl_data.py --download_species --search_term all --file_types fasta gtf --division primates --release 113 --processes 30
python download_ensembl_data.py --download_species --search_term all --file_types fasta gtf --division plants --release 60 --processes 30

# 2. Download NCBI genome files
python download_ncbi_data.py --search_term 'all' --download_species -m 1000000 --processes 20

# 3. Generate YAML file based on the genomes directory to list existing files per species as input to nf-core/references
python generate_configs.py --ensembl-dir ensembl_genomes --ncbi-dir ncbi_genomes --output-file references.yml

# 4. Run nf-core/references workflow to generate star index and intervals for species that have both fasta and gtf
nextflow run nf-core/references -profile kaust --input references.yml --outdir results/ -resume -r 0.1 --tools star,faidx,sizes

#5 . Generate genomes.config file to list the files for each species for running other nf-core workflows, similar to configs/genomes.config



