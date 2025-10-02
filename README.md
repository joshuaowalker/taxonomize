# taxonomize.py

A fast, simple tool for taxonomic classification of DNA sequences from FASTQ/FASTA files using BLAST or vsearch.

## Use Cases

- **Summarizing raw sequencing data**: Quickly identify the taxonomic composition of .fastq files to understand what's in your sequencing run
- **eDNA classification**: Classify environmental DNA reads against reference databases
- **Quality control**: Check for contamination or unexpected taxa in sequencing libraries

## Features

- Supports both BLAST and vsearch for sequence alignment
- Processes FASTA and FASTQ files
- Batch processing for memory efficiency with large read sets
- Genus or species-level reporting
- Quality filtering for FASTQ reads
- Configurable similarity thresholds and coverage requirements
- Per-read detailed output in TSV format
- Multi-threaded for performance

## Installation

### Install from GitHub

Install directly with pip (includes all Python dependencies):
```bash
pip install git+https://github.com/joshuaowalker/taxonomize.git
```

This installs the `taxonomize` command globally.

### Manual Installation

Clone the repository and install:
```bash
git clone https://github.com/joshuaowalker/taxonomize.git
cd taxonomize
pip install .
```

Or use directly without installing:
```bash
git clone https://github.com/joshuaowalker/taxonomize.git
cd taxonomize
chmod +x taxonomize.py
./taxonomize.py --help
```

### External Tools (Required)

Install at least one of:

**BLAST+** (default):
```bash
# macOS with Homebrew
brew install blast

# Ubuntu/Debian
apt-get install ncbi-blast+

# Conda
conda install -c bioconda blast
```

**vsearch** (faster alternative):
```bash
# macOS with Homebrew
brew install vsearch

# Ubuntu/Debian
apt-get install vsearch

# Conda
conda install -c bioconda vsearch
```

## Quick Start

```bash
# Basic usage with BLAST
taxonomize references.fasta reads.fastq

# Using vsearch (faster)
taxonomize references.fasta reads.fastq --search-method vsearch

# Multiple read files
taxonomize references.fasta sample1.fastq sample2.fastq sample3.fastq

# Genus-level summary with quality filtering
taxonomize references.fasta reads.fastq --level genus --min-quality 30
```

## Reference Database Format

Reference sequences must have taxonomy annotations in the header using one of these formats:

**Format 1: name= format** (recommended):
```
>seq1 name="Russula olympiana"
ATCGATCG...
>seq2 name="Russula sp. 'IN40'"
ATCGATCG...
```

**Format 2: vsearch sintax format**:
```
>seq1 g:Russula
ATCGATCG...
```

The tool extracts:
- **Genus**: First whitespace-delimited token from the name (e.g., "Russula")
- **Species**: Full name/identifier for species-level reporting (e.g., "Russula olympiana" or "Russula sp. 'IN40'")

## Usage

```
taxonomize [options] references reads [reads ...]

Positional arguments:
  references              Reference sequences file (FASTA or FASTQ)
  reads                   One or more input reads files (FASTA or FASTQ)

Search options:
  --search-method {blast,vsearch}
                          Search method to use (default: blast)
  --threads N             Number of threads (default: all available cores)
  --batch-size N          Number of reads per batch (default: 1000)

Filtering options:
  --min-identity FLOAT    Minimum sequence identity 0.0-1.0 (default: 0.95)
  --min-coverage FLOAT    Minimum query coverage 0.0-1.0 (default: 0.5)
  --min-length INT        Minimum alignment length in bp (default: 100)
  --min-quality FLOAT     Minimum mean phred quality score (default: none)
  --max-seq-length INT    Maximum sequence length in bp (default: 2000)
  --min-reads INT         Minimum reads per taxon to report (default: 2)

Output options:
  --level {genus,species} Taxonomic level for reporting (default: species)
  --max-taxa N            Maximum taxa to display per file (default: unlimited)
  --no-detailed-tsv       Don't write per-read TSV files
  --quiet                 Suppress all output except results
  --debug                 Enable debug logging
```

## Output

### Summary Output (stdout)

```
sample1.fastq	10000
  45.2%	4520	Genus species 'strain1'
  32.1%	3210	Genus species 'strain2'
  15.3%	1530	Othergenus differentspecies
  5.4%	540	Unknown
```

Format: `filename`, `total_reads`, then for each taxon: `percentage`, `count`, `taxon_name`

### Detailed Per-Read TSV

When enabled (default), creates `<filename>.taxonomize.tsv` files:

```
read_id	taxon	pident	align_length	coverage	bitscore
read1	Genus species	98.50	250	95.00	450.00
read2	Othergenus sp.	96.20	180	85.50	320.00
```

## Examples

### Quick summary of sequencing run
```bash
taxonomize database.fasta run1.fastq --max-taxa 10 --quiet
```

### High-stringency eDNA classification
```bash
taxonomize edna_refs.fasta sample.fastq \
  --min-identity 0.98 \
  --min-coverage 0.8 \
  --min-length 200 \
  --min-reads 5 \
  --level species
```

### Fast processing with vsearch
```bash
taxonomize refs.fasta reads.fastq \
  --search-method vsearch \
  --threads -1 \
  --batch-size 5000
```

### Quality-filtered genus-level report
```bash
taxonomize refs.fasta reads.fastq \
  --level genus \
  --min-quality 25 \
  --max-seq-length 1500
```

### Process multiple samples
```bash
taxonomize refs.fasta sample*.fastq --min-reads 10 --max-taxa 20
```

## Performance Tips

1. **Use vsearch for large datasets** - Generally faster than BLAST
2. **Adjust batch size** - Larger batches (5000-10000) can be faster but use more memory
3. **Use all threads** - Set `--threads -1` to use all available CPU cores
4. **Pre-filter reads** - Use quality and length filters to reduce search space

## Algorithm

1. Parse reference sequences and extract taxonomy from headers
2. Build BLAST or vsearch database from references
3. Load reads and apply quality/length filters
4. Process reads in batches:
   - Write batch to temporary FASTA
   - Search against database
   - Keep best hit per read (highest bitscore/identity)
5. Aggregate results at genus or species level
6. Filter taxa below minimum read threshold
7. Output summary and detailed TSV

## License

BSD 2-Clause License - See LICENSE file for details

Copyright (c) 2024-2025 Josh Walker

## Contributing

This is a simple, focused tool. When contributing:
- Maintain the simple wrapper design philosophy
- Ensure compatibility with standard BLAST/vsearch output formats
- Add tests for new taxonomy parsing formats
- Keep dependencies minimal
