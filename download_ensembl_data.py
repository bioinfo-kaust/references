import sys
import os
import ftplib
import logging
import argparse
from pathlib import Path
from urllib.parse import urljoin
from multiprocessing import Pool


def setup_logging():
    """Configure logging settings"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def create_ftp_connection(server):
    """Establish FTP connection"""
    try:
        ftp = ftplib.FTP(server)
        ftp.login()
        return ftp
    except Exception as e:
        logging.error(f"Failed to connect to {server}: {str(e)}")
        return None

def get_ftp_files(ftp, division, release, species, file_type, extra_dir=''):
    """Generate the appropriate FTP path based on division"""
    ftp_path = ''
    species_path = ''
    if species:
        species_path = f'/{species}{extra_dir}'
    if release!='current':
        release = f'release-{release}'
    if division == "primates":
        ftp_path = f"/pub/{release}/{file_type}{species_path}"
    else:
        ftp_path = f"/pub/{release}/{division}/{file_type}{species_path}"
    
    try:
        ftp.cwd(ftp_path)
    except Exception as e:
        logging.error(f"Error accessing {file_type} directory: {str(e)}; path tried: {ftp_path}")
    
    return ftp.nlst(), ftp_path


def download_file(ftp, remote_path, filename, local_dir):
    """Download a file from FTP server"""
    local_path = os.path.join(local_dir, filename)
    
    try:
        with open(local_path, 'wb') as f:
            ftp.retrbinary(f'RETR {remote_path}/{filename}', f.write)
        logging.info(f"Successfully downloaded: {filename}")
        return True
    except Exception as e:
        logging.error(f"Failed to download {filename}: {str(e)}")
        return False
    
def fetch_ensembl_species(ftp, species_list_file, search_term, division, release):
    """
    Fetch list of available species from Ensembl FTP server based on search criteria
    
    Args:
        ftp (ftplib.FTP): Active FTP connection object
        species_list_file (str): Path to output file where species list will be saved
        search_term (str): Term to filter species names ('all' to get all species)
        release (str, optional): Ensembl release number. 
    
    Returns:
        list: List of species names matching the search criteria or all
    """
    species_list, ftp_path = get_ftp_files(ftp, division, release, '', 'fasta', '/dna/')
    #filter by given species
    if search_term!='all':
        species_list = [x for x in species_list if search_term.lower() in x.lower()]
    
    #write to an output file
    if species_list:
        logging.info(f"Retrieved {len(species_list)} assemblies.")
        with open(species_list_file, 'w') as outfile:
            outfile.write('\n'.join(species_list))
    else:
        logging.error(f'Failed to find any species using search term {search_term}')
    # Change back to root FTP directory
    return species_list

def get_genome(species, ftp, file_types, division, release, force_replace, dna_file_ext, output_dir):
    
    release_dir = release
    file_types = file_types.split(',')
    if release!='current':
        release_dir = 'release-' + release
    output_dir = Path(output_dir) / division / species / release_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    if 'fasta' in file_types: 
        file_list, ftp_path = get_ftp_files(ftp, division, release, species.lower(), 'fasta', '/dna/')
        # Download DNA sequence file
        dna_file = f"{species}.*.{dna_file_ext}"
        matches = [f for f in file_list if f.endswith(dna_file_ext)]
        if matches:
            if not os.path.isfile(output_dir / matches[0]) or force_replace:
                logging.info(f"Downloading {ftp_path}{matches[0]}")
                download_file(ftp, ftp_path, matches[0], output_dir)
            else:
                logging.info(f"File already exists {output_dir}/{matches[0]}")
        else:
            logging.warning("DNA sequence file not found")
    
    if 'gtf' in file_types:
        file_list, ftp_path = get_ftp_files(ftp, division, release, species.lower(), 'gtf')
        matches = [f for f in file_list if f.endswith('.gtf.gz')]
        if matches:
            if not os.path.isfile(output_dir / matches[0]) or force_replace:
                logging.info(f"Downloading {ftp_path}/{matches[0]}")
                download_file(ftp, ftp_path, matches[0], output_dir)
            else:
                logging.info(f"File already exists {output_dir}/{matches[0]}")
        else:
            logging.warning("GTF file not found")
        
def main():
    parser = argparse.ArgumentParser(description='Download files from Ensembl FTP servers')
    parser.add_argument('-r', '--release', default='113', required=False, help='Ensembl release number')
    
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-s", "--species", nargs="+",
                        help="List of species names or taxon IDs")
    group.add_argument("-f", "--file", type=str,
                        help="File containing species names or taxon IDs (one per line)")
    group.add_argument("-st", "--search_term", type=str,
                        help="search term to find a species and download its files")
    parser.add_argument("-d", "--download_species", action='store_true', default=False,
                        help="Download species available in NCBI")
    parser.add_argument("-fr", "--force_replace", action='store_true', default=False,
                        help="Download and replace existing files even if it exists, especially useful to update current release")
    parser.add_argument('-fe', '--dna_file_ext', default='dna_sm.toplevel.fa.gz', help='file extension for dna fasta files default dna_sm.toplevel.fa.gz')
    parser.add_argument('-div', '--division', choices=['primates', 'plants', 'fungi', 'bacteria', 'protists', 'metazoa'],
                      default='primates', help='Ensembl division')
    parser.add_argument('-ft', '--file_types', choices=['fasta', 'gtf'], nargs="+",
                      default=['fasta', 'gtf'], help='File types to download default: fasta gtf')
    parser.add_argument("--species_list_file", default="ensembl_species_list.txt", 
                        help="File name to save the list of species names.")
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="Number of parallel downloads (default: 1)")
    parser.add_argument('-o', '--output', default='ensembl_genomes', help='Output directory')
    
    args = parser.parse_args()

    setup_logging()
    
    # Determine the appropriate FTP server
    server = "ftp.ensembl.org" if args.division == "primates" else "ftp.ensemblgenomes.org"
    
    # Connect to FTP server
    ftp = create_ftp_connection(server)
    if not ftp:
        return None
    
    # Get list of species either from comm and line or file
    species_list = []
    if args.species:
        species_list = args.species
    elif args.file and os.path.isfile(args.file):
        species_list = [line.strip() for line in open(args.file, 'r').readlines() if line.strip() if line!=""]
    elif args.search_term:
        species_list = fetch_ensembl_species(ftp, args.species_list_file, args.search_term, division=args.division, release=args.release)
    elif args.download_species:
        species_list = fetch_ensembl_species(ftp, args.species_list_file, search_term='all', division=args.division, release=args.release)
    else:
        logging.error("Please provide either \na) a list of species\nb) a file containing species names or\n c) a search term.")
    if not species_list:
        return []
    if args.search_term=='all' and not args.download_species:
        sys.exit(0)
    # Create output directory
    output_dir = Path(args.output)
    # Process species in parallel
    results = []
    if args.processes>1:
        logging.info(f'Downloading {len(species_list)} species using {args.processes} threads.')
        with Pool(processes=args.processes) as pool:
            results = pool.starmap(get_genome, [(species, ftp, ','.join(args.file_types), args.division, args.release, 
                                                 args.force_replace, args.dna_file_ext, args.output) for species in species_list])
    else:
        for species in species_list:
            results.append(get_genome(species, ftp, args.file_types, args.division, 
                                      args.release, args.force_replace, args.dna_file_ext, args.output))
    
    ftp.quit()
    return results

if __name__ == "__main__":
    main()