import csv
from sys import argv

# Define phage-related terms (case-insensitive)
phage_terms = {
    "phage", "bacteriophage", "caudoviricetes", "steitzviridae", "microviridae",
    "microvirus", "myoviridae", "prokaryotic", "siphoviridae", "siphovirus",
    "podoviridae", "podovirus", "inoviridae", "leviviridae", "herelleviridae",
    "ackermannviridae","chimeric", "cystoviridae"
}
def contains_phage(row):
    """Check if any taxonomic field contains phage-related terms."""
    fields_to_check = [
        row["Virus_Class"].lower(),
        row["Virus_Family"].lower(),
        row["Virus_Genus"].lower(),
        row["Virus_Species"].lower(),
        row["Blastn_Title"].lower(),
        row["Blastn_Species"].lower()

    ]
    return any(any(term in field for term in phage_terms) for field in fields_to_check)

reagent_terms = {"reagent", "jduan"}
def contains_reagent_terms(row):
    # Define excluded terms for Reagent
    fields_to_check = [
        row["Virus_Species"].lower(),
        row["Blastn_Title"].lower(),
        row["Blastn_Species"].lower()
    ]

    """Check if Virus_Species field contains excluded terms."""
    return any(any(term in field for term in reagent_terms) for field in fields_to_check)

def is_not_biologically_reasonable(row):
    """检查病毒宿主兼容性"""
    excluded_hosts = [
        "whale", "fish", "coral", "shrimp", "plant", "tomato", 
        "algae", "phytoplankton", "salmon"
    ]
    fields = [
        row.get("Virus_Species", "").lower(),
        row.get("Virus_Genus", "").lower(),
        row.get("Blastn_Title", "").lower(),
        row.get("Blastn_Species", "").lower()
    ]
    return any(any(term in field for term in excluded_hosts) for field in fields)



def meets_thresholds(row):
    # First check if we have any contig support at all
    has_contig_support = False
    try:
        # More robust check for contig presence
        has_contig_support = (
            row.get("Number_of_contigs", "0").strip() not in ("", "0") or
            row.get("RPM_mapping_contigs", "0").strip() not in ("", "0") or
            row.get("Mapped_Reads", "0").strip() not in ("", "0")
        )
    except (ValueError, AttributeError):
        has_contig_support = False

    # Process contig-supported cases
    if has_contig_support:
        try:
            # Get contig metrics
            rpm_contig = float(row.get("RPM_mapping_contigs", 0))
            length = float(row.get("Length", 0))
            coverage = float(row.get("Coverage", 0))
            
            read_support = (
                float(row.get("Number_of_supporting_reads", 0)) >= 3 and
                float(row.get("RPM", 0)) >= 1 and
                float(row.get("Average_identity", 0)) >= 90 and
                float(row["Average_alignment_length"]) >= 70
            )
            
            # Main contig threshold check
            contig_passes = (
                length >= 300 and 
                coverage >= 50 
            )
            
            # Either meets normal contig thresholds OR has strong read support
            return contig_passes and (rpm_contig >= 1 or read_support)
            
        except (ValueError, TypeError, AttributeError):
            return False
    
    # Process read-only cases
    else:
        try:
            return (
                float(row["Number_of_supporting_reads"]) >= 3 and
                float(row["Average_identity"]) >= 90 and
                float(row["Average_alignment_length"]) >= 70 and
                float(row["RPM"]) >= 1
            )
        except (ValueError, KeyError):
            return False

# Read the summary file
virus_summary = argv[1]
with open(virus_summary, "r") as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)
    fieldnames = reader.fieldnames

# Filter entries
retained, filtered = [], []
for row in rows:
# Always filter out phage/reagent entries, regardless of thresholds
    if contains_phage(row) or contains_reagent_terms(row) or is_not_biologically_reasonable(row):
        filtered.append(row)
    # Retain only if ALL thresholds are met
    elif meets_thresholds(row):
        retained.append(row)
    else:
        filtered.append(row)

# Write output files
virus_summary_filtered = argv[2]
with open(virus_summary_filtered, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(retained)
