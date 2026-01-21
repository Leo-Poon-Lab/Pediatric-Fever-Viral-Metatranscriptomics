import os
import glob
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO
import random

def run_blast(input_dir, output_dir, database, task="blastn", evalue=1e-5, threads=4):
    """Run BLASTn on all FASTA files in input directory"""
    # Create output directory if needed
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log_file = os.path.join(output_dir, "subset_log.txt")

    # Get all FASTA files
    fasta_files = glob.glob(os.path.join(input_dir, "*.fa"))

    # BLAST output header
    header = (
        "Query_ID\tQuery_Length\tSubject_ID\tSubject_Length\t"
        "Percentage_Identity\tNum_Identical_Matches\tPercentage_Positive\t"
        "Num_Positive_Matches\tAlignment_Length\tQuery_Coverage_Percent\t"
        "Query_Coverage_Percent_HSP\tQuery_Coverage_Percent_Unified\t"
        "Num_Mismatches\tNum_Gap_Openings\tQuery_Start\tQuery_End\t"
        "Subject_Start\tSubject_End\tNum_Gaps\tE_value\tBit_Score\t"
        "Subject_Strand\tSubject_Taxonomy_ID\tSubject_Title"
    )

    for fasta in fasta_files:
        # Generate output filename
        base_name = os.path.basename(fasta).rsplit('.', 1)[0]
        output_file = os.path.join(output_dir, f"{base_name}_blast.tsv")
        temp_output = os.path.join(output_dir, "temp_blast.tsv")

        # Read sequences
        sequences = list(SeqIO.parse(fasta, "fasta"))
        num_sequences = len(sequences)

        temp_fasta = os.path.join(output_dir, f"{base_name}_temp.fa")
        SeqIO.write(sequences, temp_fasta, "fasta")

        # Build BLAST command
        blast_cmd = [
            "blastn",
            "-query", temp_fasta,
            "-db", database,
            "-task", task,
            "-outfmt", "6 qseqid qlen sseqid slen pident nident ppos positive "
                      "length qcovs qcovhsp qcovus mismatch gapopen qstart qend "
                      "sstart send gaps evalue bitscore sstrand staxid stitle",
            "-max_target_seqs", "1",
            "-max_hsps", "1",  # Add this line
            "-evalue", str(evalue),
            "-num_threads", str(threads),
            "-out", temp_output
        ]

        # Run BLAST
        try:
            subprocess.run(blast_cmd, check=True)

            # Add header to results
            with open(temp_output, 'r') as f_in, open(output_file, 'w') as f_out:
                f_out.write(header + "\n")
                f_out.write(f_in.read())

            # Remove temporary files
            os.remove(temp_output)
            os.remove(temp_fasta)
            print(f"Successfully processed: {os.path.basename(fasta)}")

        except subprocess.CalledProcessError as e:
            print(f"Error processing {fasta}: {str(e)}")
        except Exception as e:
            print(f"Unexpected error with {fasta}: {str(e)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run BLASTn on FASTA files')
    parser.add_argument('--input_dir', required=True, help='Input directory with FASTA files')
    parser.add_argument('--output_dir', required=True, help='Output directory for BLAST results')
    parser.add_argument('--database', required=True, help='Path to BLAST database')
    parser.add_argument('--task', default='blastn', help='BLAST task type (default: blastn)')
    parser.add_argument('--evalue', type=float, default=1e-5, help='E-value threshold (default: 1e-5)')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    args = parser.parse_args()

    run_blast(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        database=args.database,
        task=args.task,
        evalue=args.evalue,
        threads=args.threads,
    )