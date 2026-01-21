import os
import subprocess
import argparse
import tempfile
from pathlib import Path

def parse_taxon_sources(taxon_ids, kraken_report_path, kaiju_output_path):
    """Determine which taxon IDs come from Kraken vs Kaiju"""
    kraken_ids = set()
    kaiju_ids = set()
    
    # Get Kraken taxon IDs from report
    if os.path.exists(kraken_report_path):
        with open(kraken_report_path) as f:
            for line in f:
                if line.strip() and not line.startswith('%'):
                    parts = line.strip().split('\t')
                    if len(parts) > 6:
                        kraken_ids.add(parts[6].strip())

    # Get Kaiju taxon IDs from kaiju.out
    if os.path.exists(kaiju_output_path):
        with open(kaiju_output_path) as f:
            for line in f:
                if line.startswith('C'):
                    parts = line.strip().split('\t')
                    if len(parts) > 2:
                        kaiju_ids.add(parts[2].strip())
    return (
        [tid for tid in taxon_ids if tid in kraken_ids],
        [tid for tid in taxon_ids if tid in kaiju_ids]
    )

def extract_reads(args, kraken_tids, kaiju_tids, output_prefix):
    """Handle paired-end read extraction with proper tool selection"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create output directory if not exists
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        
        # Process Kraken reads
        kraken_files = []
        if kraken_tids:
            # Create temporary output files
            kraken_temp_r1 = os.path.join(tmpdir, "kraken_temp_R1.fa")
            kraken_temp_r2 = os.path.join(tmpdir, "kraken_temp_R2.fa")
            
            kraken_cmd = [
                "python", "scripts/tools/KrakenTools-1.2/extract_kraken_reads.py",
                "-k", args.kraken_output,
                "-s1", args.input_r1,
                "-s2", args.input_r2,
                "-o", kraken_temp_r1,
                "-o2", kraken_temp_r2,
                "--include-children",
                "-r", args.kraken_report
                                        ]
            # Handle multiple taxids by looping through them
            for tid in kraken_tids:
                current_cmd = kraken_cmd + ["-t", tid]
                subprocess.run(current_cmd, check=True)

            # Move final outputs
            final_kraken_r1 = f"{output_prefix}_kraken_R1.fa"
            final_kraken_r2 = f"{output_prefix}_kraken_R2.fa"
            os.rename(kraken_temp_r1, final_kraken_r1)
            os.rename(kraken_temp_r2, final_kraken_r2)
            kraken_files = [final_kraken_r1, final_kraken_r2]

        # Process Kaiju reads
        kaiju_files = []
        if kaiju_tids:
            # Create read ID list for Kaiju
            kaiju_reads_file = os.path.join(tmpdir, "kaiju_reads.txt")
            with open(kaiju_reads_file, "w") as f:
                for tid in kaiju_tids:
                    awk_command = f"$1 == \"C\" && $3 == {tid} {{print $2}}"
                    subprocess.run([
                        "awk",
                        "-v", "FS=\t",  # Correct field separator specification
                        awk_command,
                        args.kaiju_output
                    ], stdout=f, check=True)

            # Extract reads using seqtk
            for end, suffix in zip([args.input_r1, args.input_r2], ["R1", "R2"]):
                output_file = f"{output_prefix}_kaiju_{suffix}.fa"
                with open(output_file, "w") as outfile:
                    # Run seqtk subseq command
                    p1 = subprocess.Popen([
                        "seqtk", "subseq",
                        end,
                        kaiju_reads_file
                    ], stdout=subprocess.PIPE, text=True)
                    
        
                    # Run seqtk seq -A command converting fq to fa
                    p2 = subprocess.Popen([
                        "seqtk", "seq", "-A"
                    ], stdin=p1.stdout, stdout=outfile, text=True)
        
                    # Close the stdout pipe of the first process
                    p1.stdout.close()
                    # Wait for both processes to complete
                    p2.communicate()
                    # Wait for both processes to complete
                    # Check for errors
                    if p1.wait() != 0:
                        raise subprocess.CalledProcessError(p1.returncode, p1.args)
                    if p2.wait() != 0:
                        raise subprocess.CalledProcessError(p2.returncode, p2.args)

                kaiju_files.append(output_file)
            

        # Merge and deduplicate results
        final_files = []
        for suffix in ["R1", "R2"]:
            # Combine Kraken and Kaiju results
            combined = f"{output_prefix}_{suffix}.fa"
            to_combine = [f for f in kraken_files + kaiju_files if suffix in f]
   
            if to_combine:
                subprocess.run(f"cat {' '.join(to_combine)} > {combined}", shell=True, text=True)
                
                # Deduplicate if needed
                final_output = f"{output_prefix}_{suffix}.combined.fa"
                if args.deduplicate and os.path.getsize(combined) > 0:
                    result = subprocess.run([
                        "seqkit", "rmdup", "-s",
                        "-o", final_output,
                        combined
                    ], text=True, capture_output=True)
                    if result.returncode != 0:
                        print(f"Error running seqkit rmdup: {result.stderr}", file=sys.stderr)
                    os.remove(combined)
                else:
                    os.rename(combined, final_output)
                
                final_files.append(final_output)
            else:
                # Create empty files if no reads
                open(f"{output_prefix}_{suffix}.combined.fa", "w").close()

        return final_files

def main():
    parser = argparse.ArgumentParser(description="Extract viral reads from Kraken/Kaiju results")
    parser.add_argument("--combined-report", required=True)
    parser.add_argument("--kraken-report", required=True)
    parser.add_argument("--kaiju-output", required=True)
    parser.add_argument("--kraken-output", required=True)
    parser.add_argument("--input-r1", required=True)
    parser.add_argument("--input-r2", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--deduplicate", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.combined_report) as f:
        headers = f.readline().strip().split()
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            record = dict(zip(headers, parts))
            
            taxon_ids = [tid for tid in record['Taxon_ID'].split(',') if tid]
            kraken_reads = int(record['Kraken_Reads'])
            kaiju_reads = int(record['Kaiju_Reads'])

            kraken_tids, kaiju_tids = parse_taxon_sources(
                taxon_ids,
                args.kraken_report,
                args.kaiju_output
            )
            # Create sanitized filename
            sample_name = record['Sample']
            taxon_name = ''.join(
                c if c.isalnum() else '_' for c in record['Taxon_Name']
            )
            output_prefix = os.path.join(
                args.output_dir,
                f"{sample_name}_{taxon_name}"
            )

            # Determine extraction strategy
            if kaiju_reads == 0 and kraken_reads > 0:
                extract_reads(args, kraken_tids, [], output_prefix)
            elif kraken_reads == 0 and kaiju_reads > 0:
                extract_reads(args, [], kaiju_tids, output_prefix)
            else:
                extract_reads(args, kraken_tids, kaiju_tids, output_prefix)

if __name__ == "__main__":
    main()