#!/usr/bin/env python3

import sys
import os
import json
import time
import random
import zipfile
import argparse
import subprocess
from pathlib import Path
import concurrent.futures
from multiprocessing import Pool

def download_genome(species, output_dir, max_retries=3, file_types='genome,gtf'):
    """
    Download genome and GTF files for a given species using NCBI datasets tool
    
    Args:
        species (str): Species name or taxon ID
        output_dir (Path): Directory to save downloaded files
        file_types: list of file types to download e.g. genome,gtf,etc
        
    """
    # Add random delay before starting download
    delay = random.uniform(1, 5)
    time.sleep(delay)
    print(f"\nProcessing {species} (after {delay:.1f}s delay)...")
    
    # Create species-specific directory
    species_dir = output_dir / species.replace(" ", "_")
    species_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if the file(s) already exist
    dataset_catalog = species_dir / "ncbi_dataset" "/data/dataset_catalog.json"
    if (dataset_catalog).exists():
        with open(dataset_catalog, 'r') as json_file:
            species_accession = json.load(json_file)
            files_index = None
            for index, element in enumerate(species_accession['assemblies']):
                try:
                    list(element.keys()).index('accession')
                    files_index = index
                    break
                except ValueError:
                    continue
            
            if files_index is None:
                print('No Accession is found in:', files_index, species_accession['assemblies'], 'Trying to download again.')  
            else:
                fasta_file = ''
                gtf_file = ''
                fasta = False
                gtf = False
                
                for file in species_accession['assemblies'][files_index]['files']:
                    file_type = file['fileType']
                    file_path = file['filePath']
                    if file_type == 'GENOMIC_NUCLEOTIDE_FASTA':
                        fasta_file = species_dir / "ncbi_dataset" "/data" / file_path
                    elif file_type == 'GTF':
                        gtf_file = species_dir / "ncbi_dataset" "/data" / file_path
                
                if fasta_file:
                    if fasta_file.exists():
                        fasta = True
                
                if gtf_file:
                    if gtf_file.exists():
                        gtf = True
                else: #indicating no annotation is available for this species (GTF is missin in the dataset report)
                    print(f"No GTF files is avilable for: {species}")
                    gtf = True
                
                if fasta and gtf:
                    print(f'Files already exist for species: {species}\n{fasta_file} {gtf_file}')
                    return (species, True)
                elif fasta and 'gtf' not in file_types:
                    print(f'Files already exist for species: {species}\n{fasta_file}')
                    return (species, True)
                elif gtf and 'genome' not in file_types:
                    print(f'Files already exist for species: {species}\n{gtf_file}')
                    return (species, True)
                else:
                    print("Files don't exist. Trying to download")
    
    #Download the files
    try:
        for attempt in range(max_retries):
            print(f'Attempt: {attempt}. Downloading: {species} to {species_dir}')
            # Download genome data using datasets command
            cmd = f"datasets download genome taxon \"{species}\" --include {file_types} --reference"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=species_dir)
            if result.returncode == 0:
                #Unzip downloaded folder that contains the files
                zip_path =  species_dir / "ncbi_dataset.zip"
                extract_path = species_dir 
                try:
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_path)
                        # Check if download was successful by looking for the data report
                        if dataset_catalog.exists():
                            # Clean up zip file
                            zip_file.unlink()
                            print(f"Successfully downloaded data for {species}")
                            zip_path.unlink()
                            return (species, True)
                except zipfile.BadZipFile as e:
                    print(f'Error in unzipping {species}. Error: {e}')
                except FileNotFoundError as e:
                    print(f'Error in unzipping {species}. Error: {e}')
            else:
                print(f"Error getting the results (code: {result.returncode}): {cmd}")
    except Exception as e:
        print(f"Error processing {species}")
    #remove directory if empty, no download is made
    try:
        species_dir.rmdir()
    except OSError:
        pass
    return (species, False)

def get_species_names(search_term, max_results=100, accepted_divisions=['all']):
    """
    Fetches a list of species names from NCBI datasets tool based on a search term.
    
    Args:
        search_term (str): The term to search for in the NCBI datasets tool.
        max_results (int): Maximum number of results to fetch.
        
    Returns:
        list: A list of species names.
    """
    try:
        cmd = f"datasets summary genome taxon all --assembly-source all --limit {max_results}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error fetching species names: {result.stderr}")
            return []
        
        data = json.loads(result.stdout)
        species_names = []
        for genome in data['assemblies']:
            if genome['organism']['name'] not in species_names:
                species_names.append(genome['organism']['name'])
        
        return species_names

    except Exception as e:
        print(f"An error occurred: {e}")
        return []

def fetch_ncbi_species(max_species, accepted_divisions, species_list_file):
    """
    Fetch species names from NCBI datasets tool.
    
    Args:
        max_species (int): Maximum number of species to fetch.
        accepted_divisions (list): List of accepted divisions.
        species_list_file (str): File name to store list of species from NCBI.
    
    Returns:
        list: List of species names.
    """
    try:
        cmd = f"datasets summary genome taxon all --assembly-source all --limit {max_species}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error fetching species names: {result.stderr}")
            return []
        
        data = json.loads(result.stdout)['reports']
        print(data)
        species_list = []
        with open(species_list_file, 'w') as species_list_outfile:
            for genome in data:
                print('genome data: ', genome)
                species_name = genome['assembly_info']['assembly_name']
                print(species_name)
                if species_name not in species_list:
                    species_list.append(species_name)
                    species_list_outfile.write(species_name + '\n')
                if len(species_list) >= max_species:
                    break
        
        return species_list

    except Exception as e:
        print(f"An error occurred: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Download genome and GTF files from NCBI for multiple species")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--species", nargs="+", 
                        help="List of species names or taxon IDs")
    group.add_argument("-f", "--file", type=str, 
                        help="File containing species names or taxon IDs (one per line)")
    group.add_argument("-st", "--search_term", type=str,
                        help="Search term to identify species name and if it exists")
    group.add_argument("-l", "--list_all_species", action='store_true', default=False,
                        help="List all species available in NCBI")
    group.add_argument("-ad", "--accepted_divisions", nargs="+", default=['all'],
                        help="""List of accepted divisions including: 'Bacteria', 'Vertebrates', 'Invertebrates', 'Viruses', 'Environmental samples', 
                            'Mammals', 'Synthetic and Chimeric', 'Unassigned', 'Phages', 'Plants and Fungi', 'Primates', 'Rodents'""")
    
    parser.add_argument("--species_list_file", default="species_list.txt", help="File name to save the list of species names when --list_all_species is enabled.")
    parser.add_argument("-d", "--download_species", action='store_true', default=False,
                        help="Download species available in NCBI, limited by --search_term or --accepted_divisions")
    parser.add_argument("-m", "--max_species", type=int, default=1,
                        help="Maximum number of species to download (default: 10k)")
    parser.add_argument("-a", "--max_attempts", type=int, default=3,
                        help="Maximum number of attempts to retry download (default: 3)")
    parser.add_argument("-ft", "--file_types", nargs="+", default=['genome', 'gtf'],
                        help="""List of assets to download e.g. genome gtf""")
    parser.add_argument("-e", "--email", type=str, default="bioinfo@kaust.edu.sa",
                        help="Email address for NCBI Entrez (default: bioinfo@kaust.edu.sa)")
    parser.add_argument("-o", "--output", type=str, default="genomes",
                        help="Output directory (default: genomes)")
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="Number of parallel downloads (default: 1)")
    args = parser.parse_args()
    
    #list all species
    if args.list_all_species:
        species = fetch_ncbi_species(args.max_species, args.accepted_divisions, args.species_list_file)
        print(f"Total species found: {len(species)}")
        for s in species[:10]:
            print(s)
            sys.exit(0)

    if args.search_term:
        print('Species found:\n', '\n'.join(get_species_names(search_term=args.search_term, max_results=args.max_species, accepted_divisions=args.accepted_divisions)))
        if not args.download_species:
            sys.exit(0)
    
    if args.max_species>10 and (args.search_term or args.list_all_species):
        print(f"Retrieving the list of {args.max_species} may take very long time..... It is faster to provide an input file with a long list of species. \nProcess Started... ")
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get list of species either from command line or file
    species_list = []
    if args.species:
        species_list = args.species
    elif args.file and os.path.isfile(args.file):
        species_list = [line.strip() for line in open(args.file, 'r').readlines() if line.strip() if line!=""]
    elif args.search_term and args.download_species:
        species_list = get_species_names(search_term=args.search_term, max_results=args.max_species, accepted_divisions=args.accepted_divisions)
    elif args.download_species and args.accepted_divisions!=['all']:
        species_list = fetch_ncbi_species(args.max_species, args.accepted_divisions, args.species_list_file)
    else:
        print("Please provide a list of species or a file containing species names.")
        sys.exit(1)
    
    # Process species in parallel
    successful = 0
    failed = 0
    failed_species = []
    
    results = []
    if args.processes>1:
        with Pool(processes=args.processes) as pool:
            results = pool.starmap(download_genome, [(spcies, output_dir, args.max_attempts, ','.join(args.file_types)) for spcies in species_list])
    else:
        for spcies in species_list:
            results.append(download_genome(spcies, output_dir, args.max_attempts, ','.join(args.file_types)))
    print(results)

if __name__ == "__main__":
    main()
