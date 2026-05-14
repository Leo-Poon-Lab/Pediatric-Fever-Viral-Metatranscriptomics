"""Main Snakemake workflow for the pediatric fever viral metatranscriptomics study."""

__maintainer__ = "Zhang Tao"
__email__ = "zhangtao@hku.hk"

import csv
from os.path import join
from pathlib import Path


configfile: "config/config.yaml"


BASE_DIR = Path(workflow.basedir).resolve()


def _config_path(key, default):
     return Path(config.get("io", {}).get(key, default))


def _output_path(key, default):
     return config.get("io", {}).get("output", {}).get(key, default)


SAMPLES_FILE = BASE_DIR / _config_path("samples_file", "config/samples.tsv")
DATA_DIR = BASE_DIR / _config_path("raw_data_dir", "resources")
LOG_DIR = BASE_DIR / _output_path("logs", "logs")

preprocessing_dir = _output_path("preprocessing", "00_preprocessing")
reads_taxa_kranken2_dir = "01_reads_taxa/kraken2"
reads_taxa_kaiju_dir = "01_reads_taxa/kaiju"
reads_taxa_kraken_kaiju_dir = "01_reads_taxa/kraken2_kaiju"
assembly_dir = _output_path("assembly", "02_assembly")
virus_identification_dir = _output_path("virus_id", "03_virus_identification")
abandunce_dir = _output_path("abundance", "04_abundance")
abandunce_consensus_dir = f"{abandunce_dir}_consensus"
reads_contigs_dir = "05_virus_reads_contigs_integration_analysis"
virus_annotation_dir = _output_path("annotation", "05_virus_annotation")
virus_component_analysis_dir = "05_virus_component_analysis"
virus_phylogenetic_dir = _output_path("phylogeny", "06_phylogeny")

LOG_DIR.mkdir(parents=True, exist_ok=True)


def load_samples(samples_file):
     if not samples_file.exists():
          raise FileNotFoundError(f"Samples sheet not found: {samples_file}")

     sample_records = {}
     with samples_file.open("r", newline="") as handle:
          reader = csv.DictReader(handle, delimiter="\t")
          required_columns = {"sample_id", "name"}
          missing = required_columns.difference(reader.fieldnames or [])
          if missing:
               missing_cols = ", ".join(sorted(missing))
               raise ValueError(f"Samples sheet is missing required columns: {missing_cols}")

          for row in reader:
               sample_id = row["sample_id"].strip()
               sample_name = row["name"].strip()
               if not sample_id or not sample_name:
                    raise ValueError("Every sample row must contain non-empty sample_id and name values")
               sample_records[sample_id] = {
                    "name": sample_name,
                    "fastq_r1": row.get("fastq_r1", "").strip(),
                    "fastq_r2": row.get("fastq_r2", "").strip(),
               }

     if not sample_records:
          raise ValueError(f"No samples were loaded from {samples_file}")

     return sample_records


sample_records = load_samples(SAMPLES_FILE)
names = {sample_id: record["name"] for sample_id, record in sample_records.items()}
SAMPLES = list(sample_records.keys())


def get_sample_record(sample_id):
     record = sample_records[sample_id]
     if record.get("fastq_r1") or record.get("fastq_r2"):
          return record

     refreshed_records = load_samples(SAMPLES_FILE)
     refreshed_record = refreshed_records.get(sample_id, record)
     sample_records[sample_id] = refreshed_record
     names[sample_id] = refreshed_record["name"]
     return refreshed_record


def _resolve_fastq_path(path_value):
     path = Path(path_value)
     if not path.is_absolute():
          path = BASE_DIR / path
     return str(path)


def _render_pattern(pattern, sample_id, sample_name, read_number):
     return (
          str(pattern)
          .replace("{sample}", sample_name)
          .replace("{sample_id}", sample_id)
          .replace("{read}", str(read_number))
     )


def get_file_path(sample_id, sample_name, read_number, base_dir=DATA_DIR):
     pattern_key = f"R{read_number}"
     configured_files = config.get("samples", {}).get("files", {})
     sample_file_overrides = configured_files.get(sample_id) or configured_files.get(sample_name) or {}
     explicit_config_fastq = sample_file_overrides.get(pattern_key)
     if explicit_config_fastq:
          return _resolve_fastq_path(explicit_config_fastq)

     record = get_sample_record(sample_id)
     explicit_fastq = record.get(f"fastq_r{read_number}")
     if explicit_fastq:
          return _resolve_fastq_path(explicit_fastq)

     naming_pattern = config.get("samples", {}).get("naming_pattern", {})
     configured_patterns = naming_pattern.get(pattern_key)
     if isinstance(configured_patterns, str):
          candidate_patterns = [configured_patterns]
     elif configured_patterns:
          candidate_patterns = list(configured_patterns)
     else:
          candidate_patterns = [
               "{sample}_R" + str(read_number) + ".fastq.gz",
               "{sample}_" + str(read_number) + ".fastq.gz",
               "{sample_id}_R" + str(read_number) + ".fastq.gz",
               "{sample_id}_" + str(read_number) + ".fastq.gz",
          ]

     candidates = [
          Path(base_dir) / sample_id / _render_pattern(pattern, sample_id, sample_name, read_number)
          for pattern in candidate_patterns
     ]
     for candidate in candidates:
          if candidate.exists():
               return str(candidate)

     return str(candidates[0])


all_outfiles = [
     str(Path(reads_contigs_dir) / "virus_summary_with_contigs_indexhopping.filtered.final.dnarna.tsv"),
]


rule all:
     input:
          all_outfiles


include: "modules/preprocessing.smk"
include: "modules/assembly.smk"
include: "modules/virus_identification_reads_based.smk"
include: "modules/virus_identification_contig_based.smk"
include: "modules/virus_abundance_calculation.smk"
include: "modules/virus_filter.smk"



