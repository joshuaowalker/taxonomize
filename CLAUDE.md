# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`taxonomize.py` is a bioinformatics tool for taxonomic classification of DNA sequences from FASTQ/FASTA files using BLAST or vsearch alignment against reference sequences.

## Running the Tool

Basic usage:
```bash
./taxonomize.py references.fasta reads.fastq
```

With common options:
```bash
./taxonomize.py references.fasta reads.fastq \
  --min-identity 0.95 \
  --min-coverage 0.5 \
  --level species \
  --threads -1
```

Search method selection:
```bash
# Use BLAST (default)
./taxonomize.py references.fasta reads.fastq --search-method blast

# Use vsearch (faster alternative)
./taxonomize.py references.fasta reads.fastq --search-method vsearch
```

## Architecture

The tool is organized around a single class `TaxonSummarizer` that:

1. **Database Creation**: Parses reference sequences, extracts taxonomy from headers (supports `g:GenusName` or `name="Genus species"` formats), creates simplified IDs, and builds BLAST/vsearch databases in temporary files

2. **Batch Processing**: Reads are processed in configurable batches (default 1000) to manage memory. Each batch is written to a temporary FASTA file and searched against the database

3. **Hit Selection**: For each query read, keeps only the best hit (highest bitscore for BLAST, highest pident for vsearch) that passes filters (min identity, coverage, alignment length)

4. **Taxonomy Mapping**: Results are aggregated at genus or species level based on `--level` flag. Taxa with fewer than `--min-reads` are filtered

5. **Output**: Produces summary statistics to stdout and optional per-read TSV files with alignment details

## Key Components

- `TaxonSummarizer`: Main class implementing context manager protocol for automatic cleanup of temporary files
- `RefInfo`: Dataclass storing taxonomy extracted from reference sequence headers
- `ReadSummary`: Dataclass containing classification results and filter statistics
- `extract_taxonomy()`: Parses taxonomy from two supported header formats
- `run_search()`: Main processing pipeline that handles both BLAST and vsearch
- `run_blast_batch()` / `run_vsearch_batch()`: Execute searches on read batches

## Dependencies

- Python 3.10+
- BioPython (Bio.SeqIO)
- tqdm for progress bars
- External tools: BLAST+ (`makeblastdb`, `blastn`) or vsearch

## Temporary Files

All temporary files (BLAST/vsearch databases, query batches) are created using `tempfile.mkstemp()` for security and tracked in `self.temp_files`. The context manager ensures cleanup on exit.

## Testing

Currently no automated tests. When adding tests, focus on:
- Taxonomy extraction from different header formats
- Batch processing with various read counts
- Filter thresholds (quality, length, identity, coverage)
- Both BLAST and vsearch search methods
