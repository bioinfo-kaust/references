import os
import json
import logging
import argparse
from pathlib import Path

def setup_logging():
    """Configure logging settings"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def generate_reference_input_config(ensembl_dir: str = "ensembl_genomes", 
                                  ncbi_dir: str = "ncbi_genomes",
                                  output_file: str = "reference_input.config",
                                  file_types_to_add = {'gtf': '.gtf', 'readme': 'README'}):
    """
    Generate input config file for nf-core/references based on genome directories.
    
    Args:
        ensembl_dir (str): Path to directory containing Ensembl genomes
        ncbi_dir (str): Path to directory containing NCBI genomes
        output_file (str): Output config file path
    """
    logging.info("Generating input config for nf-core/references...")
    
    config_content = []
    # Process both Ensembl and NCBI genome directories
    for genome_base in [ensembl_dir, ncbi_dir]:
        if not os.path.exists(genome_base):
            continue
        # Recursively walk through all subdirectories
        for division_dir in os.listdir(genome_base):
            # Skip if this is the base directory itself
            for species in os.listdir(genome_base+'/'+division_dir):
                fasta_files = []
                species_dir = Path(os.path.abspath(genome_base)) / division_dir / species
                print('browsing species_dir', species_dir)
                # Get the immediate parent directory name as the genome name
                source = ''
                
                # Recursively search for fasta/fa/fna files
                for ext in ["*.fa", "*.fasta", "*.fna"]:
                    fasta_files.extend(list(species_dir.rglob(ext)))
                if not fasta_files:
                    continue
                
                fasta_path = str(fasta_files[0])
                
                lines_to_add = []
                for file_type, file_expr in file_types_to_add.items():
                    files_found = list(species_dir.rglob(f"*{file_expr}*"))
                    if files_found:
                        file_path = str(files_found[0])
                        if 'README' in file_expr:
                            file_path = ','.join([str(x) for x in files_found])
                        lines_to_add.append(f'            {file_type}: "{file_path}"')    
                    else:
                        print(f'No {file_type} is found in {species_dir}')
                if 'ensembl' in genome_base:
                    source = f"{genome_base} {division_dir} {fasta_path.split('/')[-2]}".replace(' ', '_')
                elif 'ncbi' in genome_base:
                    source = f"{genome_base} {fasta_path.split('/')[-2]}".replace(' ', '_')
                
                config_content.extend([
                    f'- genome: "{species}"',
                    f'            fasta: "{fasta_path}"',
                    '\n'.join(lines_to_add),
                    f'            species: "{species}"',
                    f'            source: "{source}"',
                    ])
    
    # Write the config file
    with open(output_file, 'w') as f:
        f.write('\n'.join(config_content))
    
    logging.info(f"Generated {output_file}")

def generate_final_config(input_dir: str = "results", 
                         output_file: str = "configs/genomes.config"):
    """
    Generate final genomes config file after running nf-core/references.
    
    Args:
        input_dir (str): Directory containing nf-core/references output
        output_file (str): Path to output config file
    """
    logging.info("Generating final genomes config...")
    
    # Create configs directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    config_content = ["params {", "    genomes {"]
    
    if not os.path.exists(input_dir):
        logging.warn(f"Warning: Input directory {input_dir} does not exist")
        return
    
    for genome_dir in Path(input_dir).iterdir():
        if not genome_dir.is_dir():
            continue
            
        genome_name = genome_dir.name
        
        # Find relevant files
        fasta = next(genome_dir.glob("**/genome.fa*"), None)
        fai = next(genome_dir.glob("**/genome.fa*.fai"), None)
        gtf = next(genome_dir.glob("**/genes.gtf"), None)
        bed = next(genome_dir.glob("**/genes.bed"), None)
        star_index = next((d for d in genome_dir.glob("**/star") if d.is_dir()), None)
        
        if fasta:
            config_content.extend([
                f"        '{genome_name}' {{",
                f"            fasta   = '{fasta}'",
                f"            fai     = '{fai if fai else ''}'",
                f"            gtf     = '{gtf if gtf else ''}'",
                f"            bed     = '{bed if bed else ''}'",
                f"            star    = '{star_index if star_index else ''}'",
                "        }"
            ])
    
    config_content.extend(["    }", "}"])
    
    # Write the config file
    with open(output_file, 'w') as f:
        f.write('\n'.join(config_content))
    
    logging.info(f"Generated {output_file}")

def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(description='Generate reference configuration files')
    parser.add_argument('--ensembl-dir', default='ensembl_genomes',
                      help='Path to directory containing Ensembl genomes')
    parser.add_argument('--ncbi-dir', default='ncbi_genomes',
                      help='Path to directory containing NCBI genomes')
    parser.add_argument('--output-file', default='references.yml',
                      help='Output config file path')
    parser.add_argument('--file-types', type=json.loads, 
                      default='{"gtf": ".gtf", "readme": "README"}',
                      help='JSON string of file types to add (e.g., \'{"gtf": ".gtf", "readme": "README"}\')')
    
    args = parser.parse_args()
    setup_logging()
    
    # Generate input config for nf-core/references
    generate_reference_input_config(
        ensembl_dir=args.ensembl_dir,
        ncbi_dir=args.ncbi_dir,
        output_file=args.output_file,
        file_types_to_add=args.file_types
    )
    
    # Generate final config file
    #generate_final_config()

if __name__ == "__main__":
    main() 
