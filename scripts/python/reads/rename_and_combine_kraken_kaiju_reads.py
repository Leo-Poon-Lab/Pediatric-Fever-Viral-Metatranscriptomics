import os
import glob
import argparse
from pathlib import Path

def process_files(input_dir, output_dir):
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Find all R1 combined files
    r1_files = glob.glob(os.path.join(input_dir, "*_R1.combined.fa"))
    
    for r1_path in r1_files:
        # Get corresponding R2 file path
        r2_path = r1_path.replace("_R1.combined.fa", "_R2.combined.fa")
        
        if not os.path.exists(r2_path):
            print(f"Warning: Missing R2 file for {r1_path}")
            continue

        # Parse filename components
        filename = os.path.basename(r1_path)
        # Extract base name by removing _R1.combined.fa
        if not filename.endswith('_R1.combined.fa'):
            print(f"Unexpected filename format: {filename}")
            continue
        base_name = filename[:-len('_R1.combined.fa')]
        
        # Split into sample and virus using first underscore
        try:
            sample, virus = base_name.split('_', 1)
        except ValueError:
            print(f"Error splitting base name into sample and virus: {base_name}")
            continue
        
        
        # Create output filename
        output_filename = f"{sample}_{virus}_R12.combined.fa"
        output_path = os.path.join(output_dir, output_filename)
        
        # Process both R1 and R2 files
        with open(r1_path) as r1_file, open(r2_path) as r2_file, open(output_path, "w") as out_file:
            read_counter = 1
            
            # Iterate through both files simultaneously
            for (r1_header, r1_seq), (r2_header, r2_seq) in zip(
                read_fasta(r1_file), read_fasta(r2_file)):
                
                # Create new headers
                new_r1_header = f">{sample}_{virus}_{read_counter}_R1"
                new_r2_header = f">{sample}_{virus}_{read_counter}_R2"
                
                # Write to combined file
                out_file.write(f"{new_r1_header}\n{r1_seq}\n")
                out_file.write(f"{new_r2_header}\n{r2_seq}\n")
                
                read_counter += 1

def read_fasta(file_handler):
    """Generator to read FASTA files"""
    header = ""
    seq = []
    for line in file_handler:
        line = line.strip()
        if line.startswith(">"):
            if header:
                yield (header, "".join(seq))
            header = line[1:]
            seq = []
        else:
            seq.append(line)
    if header:
        yield (header, "".join(seq))

if __name__ == "__main__":
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description='Process and combine paired FASTA files')
    parser.add_argument('--input_dir', default='extracted_reads',
                        help='Input directory containing FASTA files (default: extracted_reads)')
    parser.add_argument('--output_dir', default='combined_results',
                        help='Output directory for processed files (default: combined_results)')
    
    args = parser.parse_args()
    
    # Run processing with command-line parameters
    process_files(args.input_dir, args.output_dir)