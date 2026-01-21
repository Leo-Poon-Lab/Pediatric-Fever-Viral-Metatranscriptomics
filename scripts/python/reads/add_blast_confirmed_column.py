import csv
import os
import re
import argparse

def sanitize_taxon_name(taxon_name):
    """Sanitize taxon names to match BLAST filenames."""
    return ''.join(c if c.isalnum() else '_' for c in taxon_name)

def process_blast_file(blast_path, expected_class_taxa, taxon_ids, evalue_threshold, identity_threshold, taxon_name):
    """Process a BLAST file to count confirmed reads."""
    if not os.path.exists(blast_path):
        return None, None  # Return None if BLAST file does not exist

    with open(blast_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        blast_hits = list(reader)

    total_reads = len(blast_hits)  # Total rows in BLAST file (excluding header)
    confirmed_reads = 0

    for hit in blast_hits:
        # Filter by alignment quality
        if (float(hit['E_value']) > evalue_threshold or 
            float(hit['Percentage_Identity']) < identity_threshold):
            continue

        # Get lineage of subject Lineage (pre-processed)
        subject_lineage = hit['Lineage'].split(';')
        subject_lineage_ids = hit['TaxIDs'].split(';')

        if "Herviviricetes" in expected_class_taxa:
            # If the taxon is Herviviricetes (herpesviruses), check for endogenous integration
            subject_title = hit['Subject_Title'].lower()  # Case-insensitive check
                # Check if title suggests endogenous human virus (e.g., "endogenous human herpesvirus 6")
            if ("endogenous" in subject_title and 
                    "virus" in subject_title and 
                    "human" in subject_title):
                confirmed_reads += 1  # Count as viral
            else:
                # Original logic: check expected taxa
                if any(taxid in subject_lineage_ids for taxid in taxon_ids):
                    confirmed_reads += 1
                elif any(taxon in subject_lineage for taxon in expected_class_taxa):
                    confirmed_reads += 1
        # Check if the taxon is Viruses and handle specially
        elif taxon_name == "Viruses":
            if "Viruses" in subject_lineage:
                confirmed_reads += 1
        
        else:
            # Original logic: check expected taxa
            if any(taxid in subject_lineage_ids for taxid in taxon_ids):
                confirmed_reads += 1
            elif any(taxon in subject_lineage for taxon in expected_class_taxa):
                confirmed_reads += 1

    return confirmed_reads, total_reads

def main(blast_dir, combined_report, output_report, evalue_threshold, identity_threshold):
    """Main function to process the combined report and add Blastn_confirmed column."""
    # Extract sample name from blast_dir (e.g., PMT4462 from /path/to/PMT4462/blast_results)
    # sample_name = os.path.basename(os.path.dirname(blast_dir))

    # Load original report
    with open(combined_report, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    # Process each row
    for row in rows:
        taxon_name = row['Taxon_Name']
        class_taxa = [tid.strip() for tid in row['Class'].split(',')]
        taxon_ids = [tid.strip() for tid in row['Taxon_ID'].split(',')]
        sample_name = row['Sample']
        # Find matching BLAST file (with lineage)
        sanitized_name = sanitize_taxon_name(taxon_name)
        blast_file = os.path.join(
            blast_dir,
            f"{sample_name}_{sanitized_name}_R12.combined_blast_with_lineage.tsv"
        )

        # Get confirmed and total reads from BLAST file
        confirmed, total = process_blast_file(blast_file, class_taxa, taxon_ids, evalue_threshold, identity_threshold, taxon_name)
        if confirmed is None or total is None:
            row['Blastn_confirmed'] = "NA/NA"  # Assign NA/NA if BLAST file does not exist
        else:
            row['Blastn_confirmed'] = f"{confirmed}/{total}"

    # Write updated report
    with open(output_report, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames + ['Blastn_confirmed'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Add Blastn_confirmed column to Kraken/Kaiju combined report.")
    parser.add_argument("--blast_dir", required=True, help="Directory containing BLAST files with lineage information.")
    parser.add_argument("--combined_report", required=True, help="Path to the combined Kraken/Kaiju report.")
    parser.add_argument("--output_report", required=True, help="Path to save the updated report with Blastn_confirmed column.")
    parser.add_argument("--evalue_threshold", type=float, default=1e-5, help="E-value threshold for BLAST hits (default: 1e-5).")
    parser.add_argument("--identity_threshold", type=float, default=70.0, help="Percentage identity threshold for BLAST hits (default: 70.0).")
    
    args = parser.parse_args()

    # Run the main function with command-line arguments
    main(args.blast_dir, args.combined_report, args.output_report, args.evalue_threshold, args.identity_threshold)