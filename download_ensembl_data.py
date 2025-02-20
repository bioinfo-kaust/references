import os
import ftplib
import logging
import argparse
from pathlib import Path
from urllib.parse import urljoin


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

def get_download_path(division, release, species):
    """Generate the appropriate FTP path based on division"""
    if division == "primates":
        return f"/pub/release-{release}/fasta/{species}/dna/"
    else:
        return f"/pub/release-{release}/{division}/fasta/{species}/dna/"
        

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
    
def fetch_ensembl_species(ftp, species_list_file, search_term, release):
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
    species_list = []
    try:
        # Navigate to the correct directory
        ftp_path = f"/pub/release-{release}/fasta/"
        ftp.cwd(ftp_path)
        
        # List directory contents
        if search_term=='all':
            species_list = ftp.nlst()    
        else:
            species_list = [x for x in ftp.nlst() if search_term.lower() in x.lower()]
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
    
    #write to an output file
    if species_list:
        with open(species_list_file, 'w') as outfile:
            outfile.write('\n'.join(species_list))
    else:
        print(f'Failed to find any species using search term {search_term}')
    # Change back to root FTP directory
    return species_list

def main():
    parser = argparse.ArgumentParser(description='Download files from Ensembl FTP servers')
    parser.add_argument('-r', '--release', default='113', required=False, help='Ensembl release number')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--species", nargs="+",
                        help="List of species names or taxon IDs")
    group.add_argument("-f", "--file", type=str,
                        help="File containing species names or taxon IDs (one per line)")
    group.add_argument("-st", "--search_term", type=str,
                        help="search term to find a species and download its files")
    
    parser.add_argument('-ft', '--dna_file_ext', default='dna_sm.toplevel.fa.gz', help='file extension for dna fasta files default dna_sm.toplevel.fa.gz')
    parser.add_argument('-d', '--division', choices=['primates', 'plants', 'fungi', 'bacteria', 'protists', 'metazoa'],
                      default='primates', help='Ensembl division')
    parser.add_argument("--species_list_file", default="ensembl_species_list.txt", 
                        help="File name to save the list of species names.")
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
        species_list = fetch_ensembl_species(ftp, args.species_list_file, args.search_term, args.release)
    else:
        print("Please provide either \na) a list of species\nb) a file containing species names or\n c) a search term.")
    print(ftp.pwd())
    if not species_list:
        return []
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for species in species_list:
        try:
            # Navigate to the correct directory
            ftp_path = get_download_path(args.division, args.release, species)
            ftp.cwd(ftp_path)
            print(f'now changed! to {ftp_path}')
            # List directory contents
            file_list = ftp.nlst()
            
            # Download DNA sequence file
            dna_file = f"{species}.*.{args.dna_file_ext}"
            dna_matches = [f for f in file_list if f.endswith(args.dna_file_ext)]
            print(f"downloading {dna_file}")
            if dna_matches:
                download_file(ftp, ftp_path, dna_matches[0], args.output)
            else:
                logging.warning("DNA sequence file not found")
            
            # Navigate to GTF directory and download GTF file
            if args.division == "ensembl":
                gtf_path = f"/pub/release-{args.release}/gtf/{species}"
            else:
                gtf_path = f"/pub/release-{args.release}/{args.division}/gtf/{species}"
            
            try:
                ftp.cwd(gtf_path)
                file_list = ftp.nlst()
                gtf_matches = [f for f in file_list if f.endswith('.gtf.gz')]
                if gtf_matches:
                    download_file(ftp, gtf_path, gtf_matches[0], output_dir)
                else:
                    logging.warning("GTF file not found")
            except Exception as e:
                logging.error(f"Error accessing GTF directory: {str(e)}")
            
        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
    
    ftp.quit()

if __name__ == "__main__":
    main()