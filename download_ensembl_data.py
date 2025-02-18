import ftplib
import os
import argparse
from urllib.parse import urljoin
import logging

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

def main():
    parser = argparse.ArgumentParser(description='Download files from Ensembl FTP servers')
    parser.add_argument('-r', '--release', default='113', required=False, help='Ensembl release number')
    parser.add_argument('-s', '--species', required=True, help='Species name (e.g., homo_sapiens)')
    parser.add_argument('-ft', '--dna_file_ext', default='dna_sm.toplevel.fa.gz', help='file extension for dna fasta files default dna_sm.toplevel.fa.gz')
    parser.add_argument('-d', '--division', choices=['primates', 'plants', 'fungi', 'bacteria', 'protists', 'metazoa'],
                      default='primates', help='Ensembl division')
    parser.add_argument('-o', '--output', default='genomes', help='Output directory')
    
    args = parser.parse_args()

    setup_logging()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Determine the appropriate FTP server
    server = "ftp.ensembl.org" if args.division == "primates" else "ftp.ensemblgenomes.org"
    
    # Connect to FTP server
    ftp = create_ftp_connection(server)
    if not ftp:
        return
    try:
        # Navigate to the correct directory
        ftp_path = get_download_path(args.division, args.release, args.species)
        print(ftp_path)
        ftp.cwd(ftp_path)
        
        # List directory contents
        file_list = ftp.nlst()
        
        # Download DNA sequence file
        dna_file = f"{args.species}.*.{args.dna_file_ext}"
        dna_matches = [f for f in file_list if f.endswith(args.dna_file_ext)]
        print(f"downloading {dna_file}")
        if dna_matches:
            download_file(ftp, ftp_path, dna_matches[0], args.output)
        else:
            logging.warning("DNA sequence file not found")
        
        # Navigate to GTF directory and download GTF file
        if args.division == "ensembl":
            gtf_path = f"/pub/release-{args.release}/gtf/{args.species}"
        else:
            gtf_path = f"/pub/release-{args.release}/{args.division}/gtf/{args.species}"
        
        try:
            ftp.cwd(gtf_path)
            file_list = ftp.nlst()
            gtf_matches = [f for f in file_list if f.endswith('.gtf.gz')]
            if gtf_matches:
                download_file(ftp, gtf_path, gtf_matches[0], args.output)
            else:
                logging.warning("GTF file not found")
        except Exception as e:
            logging.error(f"Error accessing GTF directory: {str(e)}")
        
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
    finally:
        ftp.quit()

if __name__ == "__main__":
    main()