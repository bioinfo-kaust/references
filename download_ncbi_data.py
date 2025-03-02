#!/usr/bin/env python3
import os
import json
import time
import random
import zipfile
import hashlib
import logging
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool

def setup_logging():
    """Configure logging settings"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def validate_checksum(file_path, expected_md5):
    """
    Validate file using MD5 checksum
    
    Args:
        file_path (Path): Path to the file
        expected_md5 (str): Expected MD5 checksum
    
    Returns:
        bool: True if checksum matches, False otherwise
    """
    if not file_path.exists():
        return False
        
    md5_hash = hashlib.md5()
    with open(file_path, 'rb') as f:
        # Read file in chunks to handle large files efficiently
        for chunk in iter(lambda: f.read(4096), b''):
            md5_hash.update(chunk)
    return md5_hash.hexdigest() == expected_md5

def get_file_paths(dataset_catalog, file_types):
    
    if not dataset_catalog.exists():
        return {}
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
    
    base_dir = dataset_catalog.parent
    files = {}
    for file in species_accession['assemblies'][files_index]['files']:
        #only check required files
        file_type = file['fileType'].lower()
        if 'fasta' not in file_type and file_type not in file_types:
            continue
        if (base_dir / file['filePath']).exists():    
            files[file_type] = file['filePath']
    
    for u_file_type in file_types:
        if u_file_type not in files.keys() and 'fasta' not in u_file_type:
            logging.warn(f'No file was found for {u_file_type}')
            continue
        

    return files
    
def download_genome(species, output_dir, max_retries=3, file_types=['genome', 'gtf']):
    """
    Download genome and GTF files for a given species using NCBI datasets tool
    
    Args:
        species (str): Species name or taxon ID
        output_dir (Path): Directory to save downloaded files
        file_types: list of file types to download e.g. genome,gtf,etc
        
    """
    # Create species-specific directory
    species_dir = output_dir / species.replace(" ", "_")
    species_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if the file(s) already exist
    dataset_catalog = species_dir / "ncbi_dataset" "/data/dataset_catalog.json"
    md5_path = species_dir / "md5sum.txt"
    files = []
    if dataset_catalog.exists():
        files = get_file_paths(dataset_catalog, file_types)
        if files:
            logging.info(f'Files already exist for {species}')
            return (species, True)
    
    #Download the files
    # Add random delay before starting download
    delay = random.uniform(1, 5)
    time.sleep(delay)
    logging.info(f"\nProcessing {species} (after {delay:.1f}s delay)...")
    
    try:
        for attempt in range(max_retries):
            logging.info(f'Attempt: {attempt}. Downloading: {species} to {species_dir}')
            # Download genome data using datasets command
            file_types_to_get = ','.join(file_types).replace('fasta', 'genome')
            cmd = f"datasets download genome taxon \"{species}\" --include {file_types_to_get} --reference"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=species_dir)
            if result.returncode == 0:
                #Unzip downloaded folder that contains the files
                zip_file =  species_dir / "ncbi_dataset.zip"
                extract_path = species_dir
                try:
                    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                        zip_ref.extractall(extract_path)
                        # Validate checksums after extraction
                        if dataset_catalog.exists():
                            md5_hashes = {'/'.join(x.strip().split(' ')[-1].split('/')[-2:]):x.strip().split(' ')[0]
                                          for x in open(md5_path, 'r').readlines()}
                            
                            files = get_file_paths(dataset_catalog, file_types)
                            for file_type, file_path in files.items():
                                files_valid = True
                                file_full_path = species_dir / "ncbi_dataset" / "data" / file_path
                                if not validate_checksum(file_full_path, md5_hashes[file_path]):
                                    logging.warn(f"Checksum validation failed for {file_path}")
                                    files_valid = False
                                    break
                            
                                if files_valid:
                                    logging.info(f"Successfully downloaded and validated data for {species}")
                                    zip_file.unlink()
                                    return (species, True)
                                else:
                                    logging.warn(f"Checksum validation failed for {species}. Retrying...")
                                    zip_file.unlink()
                                    continue
                except zipfile.BadZipFile as e:
                    logging.error(f'Error in unzipping {species}. Error: {e}')
                except FileNotFoundError as e:
                    logging.error(f'Error in unzipping {species}. Error: {e}')
            else:
                logging.error(f"Error getting the results (code: {result.returncode}): {cmd}")
    except Exception as e:
        logging.error(f"Error processing {species}", e)
    #remove directory if empty, no download is made
    try:
        species_dir.rmdir()
    except OSError:
        pass
    return (species, False)

def fetch_ncbi_species(max_species, species_list_file=None, search_term='all', assembly_source='all', reference=True):
    """
    Fetch species names from NCBI datasets tool.
    
    Args:
        max_species (int): Maximum number of species to fetch.
        species_list_file (str): File name to store list of species from NCBI.
    
    Returns:
        list: List of species names.
    """
    try:
        cmd = f"datasets summary genome taxon '{search_term}' --assembly-source {assembly_source} --limit {max_species} --as-json-lines"
        if reference: #limit to reference genomes only
            cmd+=' --reference'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"Error fetching species names: {result.stderr}")
            return []
        data_sets = result.stdout.split('\n')
        species_list = []
        
        for genome in data_sets:
            if genome:
                species_name = json.loads(genome)['organism']['organism_name'].replace("'", "").replace("[", "").replace("]", "")
                #species_name = genome['assembly_info']['assembly_name']
                if species_name not in species_list:
                    species_list.append(species_name)
                if len(species_list) >= max_species:
                    break
    except Exception as e:
        logging.error(f"An error occurred: {e}")
    
    if species_list_file and species_list:
        with open(species_list_file, 'w') as species_list_outfile:
            species_list_outfile.write('\n'.join(list(set(species_list))))
    return species_list

def main():
    parser = argparse.ArgumentParser(description="Download genome and GTF files from NCBI for multiple species")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-s", "--species", nargs="+", 
                        help="List of species names or taxon IDs")
    group.add_argument("-f", "--file", type=str, 
                        help="File containing species names or taxon IDs (one per line)")
    parser.add_argument("-st", "--search_term", type=str, default='all',
                        help="Search term to identify species name and if it exists")
    
    parser.add_argument("-ar", "--assembly_source", type=str, default='all',
                        help="Assembly resources to use (either RefSeq or GenBank or all (default))")
    parser.add_argument("-r", "--reference_only", action='store_true', default=True,
                        help="All consider reference genomes, default True")
    parser.add_argument("-m", "--max_species", type=int, default=1,
                        help="Maximum number of species to download (default: 10k)")
    
    parser.add_argument("-d", "--download_species", action='store_true', default=False,
                        help="Download species available in NCBI")
    parser.add_argument("-a", "--max_attempts", type=int, default=3,
                        help="Maximum number of attempts to retry download (default: 3)")
    parser.add_argument("-ft", "--file_types", nargs="+", default=['fasta', 'gtf'],
                        help="""List of assets to download e.g. genome gtf""")
    
    parser.add_argument("--species_list_file", default="ncbi_species_list.txt", 
                        help="File name to save the list of species names.")
    parser.add_argument("-o", "--output", type=str, default="ncbi_genomes",
                        help="Output directory (default: ncbi_genomes)")
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="Number of parallel downloads (default: 1)")
    args = parser.parse_args()
    setup_logging()
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get list of species either from command line or file
    species_list = []
    if args.species:
        species_list = args.species
    elif args.file and os.path.isfile(args.file):
        species_list = [line.strip() for line in open(args.file, 'r').readlines() if line.strip() if line!=""]
    elif args.search_term or args.download_species:
        species_list = fetch_ncbi_species(args.max_species, args.species_list_file, args.search_term, args.assembly_source, args.reference_only)
    else:
        logging.info("Please provide either \na) a list of species\nb) a file containing species names or\n c) a search term.")

    if not species_list:
        return []
    
    # Process species in parallel
    results = []
    if args.processes>1:
        with Pool(processes=args.processes) as pool:
            results = pool.starmap(download_genome, [(spcies, output_dir, args.max_attempts, ','.join(args.file_types)) for spcies in species_list])
    else:
        for spcies in species_list:
            results.append(download_genome(spcies, output_dir, args.max_attempts, args.file_types))
    return results

if __name__ == "__main__":
    main()
