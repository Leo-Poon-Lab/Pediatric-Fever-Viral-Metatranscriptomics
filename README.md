# Pediatric-Fever-Viral-Metatranscriptomics

Snakemake workflow and analysis utilities supporting the manuscript "Elucidating the Viral Aetiology of Unexplained Pediatric Febrile Illnesses via Metatranscriptomic Sequencing".

## Study Summary

This repository accompanies a hospital-based metatranscriptomic study of 490 pediatric febrile illness cases collected in Hong Kong between 2020 and 2024. The computational workflow supports unbiased viral discovery from paired-end total RNA sequencing, including both reads-based and contig-based evidence integration.

At a high level, the active workflow performs:

- read preprocessing: quality trimming, host depletion, rRNA depletion
- reads-based viral detection: Kraken2, Kaiju, BLAST-based confirmation
- contig-based viral detection: MEGAHIT assembly, Diamond BLASTx, BLASTn confirmation
- abundance and coverage summarization
- reads/contigs integration and final study-specific filtering

This repository is intended to expose the analysis logic used in the manuscript. It does not include raw sequencing data, local reference databases, or full intermediate outputs.

## What This Repository Includes

- `Snakefile`: main Snakemake entrypoint
- `modules/`: workflow rules grouped by analysis stage
- `scripts/python/`: custom Python utilities used by the workflow
- `scripts/shell/`: supporting shell helpers
- `envs/`: Conda environments for workflow execution
- `config/config.example.yaml`: portable configuration template
- `config/samples.example.tsv`: public sample sheet template
- `config/patient_metadata.example.tsv`: public metadata template
- `TOP_JOURNAL_CODE_REVIEW.md`: engineering audit notes for the revised codebase

## What This Repository Deliberately Excludes

- raw FASTQ files
- patient-level private metadata files used for local execution
- local machine configuration files
- large external databases such as BLAST, Kraken2, Kaiju, and host reference indexes
- large intermediate workflow outputs

## Public Release Structure

For a public GitHub release, keep tracked files limited to reusable code, templates, and documentation. Real runtime files such as `config/config.yaml`, `config/samples.tsv`, `config/patient_metadata.tsv`, `data/raw/`, `.snakemake/`, and large result directories should remain local only.

Detailed release recommendations are documented in `GITHUB_RELEASE_CHECKLIST.md`.

## Repository Layout

```text
.
├── Snakefile
├── README.md
├── TOP_JOURNAL_CODE_REVIEW.md
├── config/
│   ├── config.example.yaml
│   ├── patient_metadata.example.tsv
│   └── samples.example.tsv
├── envs/
├── modules/
├── scripts/
└── data/
    └── raw/                  # local only, ignored from Git
```

## Configuration

Create local runtime files from the public templates:

```bash
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/patient_metadata.example.tsv config/patient_metadata.tsv
```

Then edit:

- `config/config.yaml` for local database paths and runtime locations
- `config/samples.tsv` for the samples to process
- `config/patient_metadata.tsv` for the study metadata used by downstream filtering

## Input Layout

By default, raw reads are expected under the directory configured by `io.raw_data_dir` in `config/config.yaml`.

Example layout:

```text
data/raw/
    SAMPLE_001/
        SAMPLE_001_R1.fastq.gz
        SAMPLE_001_R2.fastq.gz
```

The workflow accepts both `*_R1.fastq.gz` / `*_R2.fastq.gz` and `*_1.fastq.gz` / `*_2.fastq.gz` naming schemes.

If FASTQ files live outside the repository, they can be referenced explicitly through optional `fastq_r1` and `fastq_r2` columns in `config/samples.tsv`, or by using the `samples.files` mapping in `config/config.yaml`.

## Setup

Create a Snakemake runtime environment:

```bash
conda env create -f envs/snakemake.yaml -n pfvm-snakemake
conda activate pfvm-snakemake
```

## Usage

Dry-run the active workflow:

```bash
snakemake -n --use-conda
```

Run the active target locally:

```bash
snakemake --use-conda --cores 16
```

Run a specific output:

```bash
snakemake --use-conda --cores 16 00_preprocessing/01.trimmed/SAMPLE_001.trimmed_1.fq
```

## Reproducibility Notes

- Active workflow rules now use real Conda environment files rather than placeholder environment names.
- The workflow supports external FASTQ paths so large raw data do not need to be duplicated into the repository.
- Metadata has been unified to a single runtime file, `config/patient_metadata.tsv`, with `PCR` used as the active downstream column.
- The current codebase is optimized for this manuscript's cohort and filtering strategy rather than for broad clinical deployment without adaptation.
