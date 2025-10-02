#!/usr/bin/env python3.10

"""
taxonomize.py - Summarize read alignments at taxonomic level from FASTQ/FASTA files

Copyright (c) 2024-2025 Josh Walker

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Processes input read files and provides taxonomic summaries (genus or species level)
using BLAST to identify best matches against reference sequences.
"""

import argparse
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import multiprocessing
from tqdm import tqdm

# Constants
MAX_TARGET_SEQS = 10  # Maximum BLAST hits to consider per query
DEFAULT_BATCH_SIZE = 1000  # Default number of reads per batch
DEFAULT_MIN_IDENTITY = 0.95  # Default minimum identity threshold
DEFAULT_MIN_COVERAGE = 0.5  # Default minimum query coverage
DEFAULT_MIN_LENGTH = 100  # Default minimum alignment length
DEFAULT_MIN_READS = 2  # Default minimum reads per taxon to report
DEFAULT_MAX_SEQ_LENGTH = 2000  # Default maximum sequence length
BLAST_EVALUE = 1e-5  # E-value threshold for BLAST searches


@dataclass
class ReadSummary:
    """Summary statistics for a set of reads."""
    total_reads: int
    taxon_counts: Counter
    unknown_count: int
    read_details: Dict[str, Tuple[str, float, int, float, float]]  # query_id -> (taxon, pident, length, coverage, bitscore)
    filtered_low_quality: int = 0  # Reads filtered by quality score
    filtered_too_long: int = 0  # Reads filtered by sequence length


@dataclass
class RefInfo:
    """Information about a reference sequence."""
    original_id: str
    genus: Optional[str]
    species: Optional[str]  # Full species name/temp code
    description: str


class TaxonSummarizer:
    """Taxonomic summarizer using BLAST for read classification.

    Supports context manager protocol for automatic cleanup.
    """

    def __init__(self, similarity_threshold: float = DEFAULT_MIN_IDENTITY,
                 max_taxa: Optional[int] = None, threads: int = 1,
                 batch_size: int = DEFAULT_BATCH_SIZE,
                 taxonomic_level: str = "species",
                 min_reads: int = DEFAULT_MIN_READS,
                 min_quality: Optional[float] = None,
                 max_seq_length: int = DEFAULT_MAX_SEQ_LENGTH,
                 search_method: str = "blast"):
        """Initialize taxonomic summarizer.

        Args:
            similarity_threshold: Minimum identity threshold (0.0-1.0)
            max_taxa: Maximum number of taxa to show, None for unlimited
            threads: Number of threads for BLAST/vsearch
            batch_size: Number of reads to process per batch
            taxonomic_level: Either "genus" or "species" for reporting level
            min_reads: Minimum reads per taxon to report
            min_quality: Minimum mean phred quality score (None to disable filtering)
            max_seq_length: Maximum sequence length to process
            search_method: Either "blast" or "vsearch" for sequence search

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        # Validate inputs
        if not 0.0 <= similarity_threshold <= 1.0:
            raise ValueError(f"similarity_threshold must be between 0.0 and 1.0, got {similarity_threshold}")
        if max_taxa is not None and max_taxa < 1:
            raise ValueError(f"max_taxa must be >= 1 or None, got {max_taxa}")
        if threads < 1:
            raise ValueError(f"threads must be >= 1, got {threads}")
        if batch_size < 1:
            raise ValueError(f"batch_size must be >= 1, got {batch_size}")
        if taxonomic_level not in ("genus", "species"):
            raise ValueError(f"taxonomic_level must be 'genus' or 'species', got {taxonomic_level}")
        if min_reads < 1:
            raise ValueError(f"min_reads must be >= 1, got {min_reads}")
        if min_quality is not None and min_quality < 0:
            raise ValueError(f"min_quality must be >= 0 or None, got {min_quality}")
        if max_seq_length < 1:
            raise ValueError(f"max_seq_length must be >= 1, got {max_seq_length}")
        if search_method not in ("blast", "vsearch"):
            raise ValueError(f"search_method must be 'blast' or 'vsearch', got {search_method}")

        self.min_identity = similarity_threshold
        self.max_taxa = max_taxa
        self.threads = threads
        self.batch_size = batch_size
        self.taxonomic_level = taxonomic_level
        self.min_reads = min_reads
        self.min_quality = min_quality
        self.max_seq_length = max_seq_length
        self.search_method = search_method

        # Reference sequence info
        self.ref_info: Dict[str, RefInfo] = {}  # simplified_id -> RefInfo
        self.blast_db: Optional[str] = None
        self.vsearch_db: Optional[Path] = None
        self.temp_files: List[Path] = []

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with automatic cleanup."""
        self.cleanup()
        return False

    def cleanup(self):
        """Remove all temporary files."""
        if self.blast_db:
            for suffix in ['.nhr', '.nin', '.nsq']:
                try:
                    Path(self.blast_db + suffix).unlink()
                except FileNotFoundError:
                    pass
            self.blast_db = None

        if self.vsearch_db:
            try:
                self.vsearch_db.unlink()
            except FileNotFoundError:
                pass
            self.vsearch_db = None

        for temp_file in self.temp_files:
            try:
                temp_file.unlink()
            except FileNotFoundError:
                pass
        self.temp_files = []

    def extract_taxonomy(self, description: str) -> Tuple[Optional[str], Optional[str]]:
        """Extract genus and species from reference sequence description.

        Supports two formats:
        1. vsearch --sintax format: g:GenusName
        2. name= format: name="Genus species 'code'" or name="Genusaceae sp. 'code'"

        Returns:
            Tuple of (genus, species) where genus is first token and species is full name
        """
        # Try name="..." format first
        name_match = re.search(r'name="([^"]+)"', description)
        if name_match:
            full_name = name_match.group(1)
            # Extract genus (first whitespace-delimited token)
            tokens = full_name.split(None, 1)  # Split on first whitespace
            genus = tokens[0] if tokens else None
            return genus, full_name

        # Try vsearch --sintax g: format
        genus_match = re.search(r'g:(\w+)', description)
        if genus_match:
            genus = genus_match.group(1)
            return genus, genus

        return None, None

    def _get_sequence_format(self, file_path: Path) -> str:
        """Determine sequence file format from extension.

        Args:
            file_path: Path to sequence file

        Returns:
            'fasta' or 'fastq'
        """
        return "fasta" if file_path.suffix.lower() in {'.fa', '.fasta'} else "fastq"

    def _calculate_mean_quality(self, record: SeqRecord) -> Optional[float]:
        """Calculate mean phred quality score for a read.

        Args:
            record: SeqRecord with quality scores

        Returns:
            Mean quality score, or None if no quality scores available
        """
        if not hasattr(record, 'letter_annotations') or 'phred_quality' not in record.letter_annotations:
            return None

        qualities = record.letter_annotations['phred_quality']
        if not qualities:
            return None

        return sum(qualities) / len(qualities)

    def _filter_by_quality(self, records: List[SeqRecord]) -> Tuple[List[SeqRecord], int]:
        """Filter reads by mean quality score if min_quality is set.

        Args:
            records: List of SeqRecord objects

        Returns:
            Tuple of (filtered_records, num_filtered)
        """
        if self.min_quality is None:
            return records, 0

        filtered = []
        num_filtered = 0

        for record in records:
            mean_qual = self._calculate_mean_quality(record)
            # If no quality scores available (FASTA), pass through
            if mean_qual is None or mean_qual >= self.min_quality:
                filtered.append(record)
            else:
                num_filtered += 1

        return filtered, num_filtered

    def _filter_by_length(self, records: List[SeqRecord]) -> Tuple[List[SeqRecord], int]:
        """Filter reads by maximum sequence length.

        Args:
            records: List of SeqRecord objects

        Returns:
            Tuple of (filtered_records, num_filtered)
        """
        filtered = []
        num_filtered = 0

        for record in records:
            if len(record.seq) <= self.max_seq_length:
                filtered.append(record)
            else:
                num_filtered += 1

        return filtered, num_filtered

    def create_blast_db(self, ref_file: Path) -> str:
        """Create a BLAST database from reference sequences.

        Args:
            ref_file: Path to reference sequences file

        Returns:
            Path to BLAST database prefix

        Raises:
            ValueError: If no valid references found or BLAST fails
        """
        format = self._get_sequence_format(ref_file)

        # Create temporary FASTA file with simplified headers (secure)
        fd, temp_fasta_path = tempfile.mkstemp(suffix='.fasta')
        temp_fasta = Path(temp_fasta_path)
        self.temp_files.append(temp_fasta)

        # Close the file descriptor from mkstemp
        import os
        os.close(fd)

        with open(temp_fasta, 'w') as f:
            for idx, record in enumerate(SeqIO.parse(str(ref_file), format), start=1):
                if len(record.seq) == 0:
                    continue

                # Create simplified ID and store mapping
                simple_id = f"ref_{idx}"
                genus, species = self.extract_taxonomy(record.description)

                if genus:
                    self.ref_info[simple_id] = RefInfo(
                        original_id=record.id,
                        genus=genus,
                        species=species,
                        description=record.description
                    )

                    # Write sequence with simplified header
                    f.write(f">{simple_id}\n{str(record.seq)}\n")

        if not self.ref_info:
            raise ValueError("No valid reference sequences with genus information found")

        # Create BLAST database from temporary FASTA (secure)
        fd_db, db_prefix = tempfile.mkstemp(prefix='blast_db_')
        os.close(fd_db)
        db_prefix_path = Path(db_prefix)
        # Remove the file created by mkstemp, makeblastdb will create its own files
        db_prefix_path.unlink()

        cmd = [
            'makeblastdb',
            '-in', str(temp_fasta),
            '-dbtype', 'nucl',
            '-out', str(db_prefix_path)
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.debug(f"makeblastdb output: {result.stdout}")
        except FileNotFoundError:
            raise RuntimeError("makeblastdb not found. Please ensure BLAST+ is installed and in PATH.")
        except subprocess.CalledProcessError as e:
            logging.error(f"makeblastdb failed: {e.stderr}")
            raise RuntimeError(f"makeblastdb failed: {e.stderr}")

        logging.info(f"Created temporary BLAST database with {len(self.ref_info)} reference sequences")
        self.blast_db = str(db_prefix_path)
        return str(db_prefix_path)

    def create_vsearch_db(self, ref_file: Path) -> Path:
        """Create a vsearch UDB database from reference sequences.

        Args:
            ref_file: Path to reference sequences file

        Returns:
            Path to vsearch UDB database

        Raises:
            ValueError: If no valid references found or vsearch fails
        """
        format = self._get_sequence_format(ref_file)

        # Create temporary FASTA file with simplified headers (secure)
        fd, temp_fasta_path = tempfile.mkstemp(suffix='.fasta')
        temp_fasta = Path(temp_fasta_path)
        self.temp_files.append(temp_fasta)

        # Close the file descriptor from mkstemp
        import os
        os.close(fd)

        with open(temp_fasta, 'w') as f:
            for idx, record in enumerate(SeqIO.parse(str(ref_file), format), start=1):
                if len(record.seq) == 0:
                    continue

                # Create simplified ID and store mapping
                simple_id = f"ref_{idx}"
                genus, species = self.extract_taxonomy(record.description)

                if genus:
                    self.ref_info[simple_id] = RefInfo(
                        original_id=record.id,
                        genus=genus,
                        species=species,
                        description=record.description
                    )

                    # Write sequence with simplified header
                    f.write(f">{simple_id}\n{str(record.seq)}\n")

        if not self.ref_info:
            raise ValueError("No valid reference sequences with genus information found")

        # Create vsearch UDB database from temporary FASTA (optional but faster)
        fd_db, db_path = tempfile.mkstemp(suffix='.udb')
        os.close(fd_db)
        db_path_obj = Path(db_path)
        db_path_obj.unlink()  # Remove file, vsearch will create it

        cmd = [
            'vsearch',
            '--makeudb_usearch', str(temp_fasta),
            '--output', str(db_path_obj)
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.debug(f"vsearch makeudb output: {result.stderr}")
            # Track UDB file for cleanup
            self.temp_files.append(db_path_obj)
        except FileNotFoundError:
            # If vsearch not found, just use the FASTA directly
            logging.info(f"vsearch not found or UDB creation failed, will use FASTA directly")
            self.vsearch_db = temp_fasta
            return temp_fasta
        except subprocess.CalledProcessError as e:
            logging.warning(f"vsearch UDB creation failed: {e.stderr}, will use FASTA directly")
            self.vsearch_db = temp_fasta
            return temp_fasta

        logging.info(f"Created temporary vsearch UDB database with {len(self.ref_info)} reference sequences")
        self.vsearch_db = db_path_obj
        return db_path_obj

    def run_blast_batch(self, batch_records: List[SeqRecord], min_coverage: float,
                        min_length: int) -> Dict[str, Tuple[str, float, int, float, float]]:
        """Run BLAST on a batch of reads and return best hits.

        Args:
            batch_records: List of SeqRecord objects to process
            min_coverage: Minimum query coverage (0.0-1.0)
            min_length: Minimum alignment length

        Returns:
            Dictionary mapping query_id -> (taxon, pident, length, coverage, bitscore)
            where taxon is genus or species depending on self.taxonomic_level

        Raises:
            RuntimeError: If BLAST is not found or fails
        """
        if not batch_records:
            return {}

        # Create temporary query file for this batch (secure)
        import os
        fd, temp_query_path = tempfile.mkstemp(suffix='.fasta')
        temp_query = Path(temp_query_path)
        self.temp_files.append(temp_query)
        os.close(fd)

        with open(temp_query, 'w') as f:
            SeqIO.write(batch_records, f, 'fasta')

        # Run BLAST on this batch
        cmd = [
            'blastn',
            '-query', str(temp_query),
            '-db', self.blast_db,
            '-outfmt', '6 qacc sacc pident length qcovs bitscore',
            '-evalue', str(BLAST_EVALUE),
            '-perc_identity', str(self.min_identity * 100),
            '-qcov_hsp_perc', str(min_coverage * 100),
            '-max_target_seqs', str(MAX_TARGET_SEQS),
            '-num_threads', str(self.threads)
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.debug(f"BLAST batch found {len(result.stdout.splitlines())} hits")
        except FileNotFoundError:
            raise RuntimeError("blastn not found. Please ensure BLAST+ is installed and in PATH.")
        except subprocess.CalledProcessError as e:
            logging.error(f"BLAST search failed: {e.stderr}")
            raise RuntimeError(f"BLAST search failed: {e.stderr}")

        # Process BLAST results for this batch
        batch_best_hits = {}

        for line in result.stdout.splitlines():
            if not line.strip():
                continue

            query_id, sseqid, pident, length, coverage, bitscore = line.split('\t')
            pident = float(pident)
            length = int(length)
            coverage = float(coverage)
            bitscore = float(bitscore)

            # Apply filters
            if (pident >= (self.min_identity * 100) and
                    length >= min_length):

                if sseqid in self.ref_info and self.ref_info[sseqid].genus:
                    ref = self.ref_info[sseqid]
                    # Use genus or species depending on taxonomic_level
                    taxon = ref.species if self.taxonomic_level == "species" else ref.genus
                    # Keep hit if it's better than what we've seen for this read
                    if (query_id not in batch_best_hits or
                            bitscore > batch_best_hits[query_id][4]):
                        batch_best_hits[query_id] = (taxon, pident, length, coverage, bitscore)
                        if logging.getLogger().isEnabledFor(logging.DEBUG):
                            logging.debug(
                                f"{query_id}: {taxon} "
                                f"(pident={pident:.1f}%, len={length}, "
                                f"cov={coverage:.1f}%, score={bitscore:.1f})"
                            )

        return batch_best_hits

    def run_vsearch_batch(self, batch_records: List[SeqRecord], min_coverage: float,
                          min_length: int) -> Dict[str, Tuple[str, float, int, float, float]]:
        """Run vsearch on a batch of reads and return best hits.

        Args:
            batch_records: List of SeqRecord objects to process
            min_coverage: Minimum query coverage (0.0-1.0)
            min_length: Minimum alignment length

        Returns:
            Dictionary mapping query_id -> (taxon, pident, length, coverage, bitscore)
            where taxon is genus or species depending on self.taxonomic_level
            Note: bitscore will always be 0.0 for vsearch (not computed)

        Raises:
            RuntimeError: If vsearch is not found or fails
        """
        if not batch_records:
            return {}

        # Create temporary query file for this batch (secure)
        import os
        fd, temp_query_path = tempfile.mkstemp(suffix='.fasta')
        temp_query = Path(temp_query_path)
        self.temp_files.append(temp_query)

        os.close(fd)

        with open(temp_query, 'w') as f:
            SeqIO.write(batch_records, f, 'fasta')

        # Run vsearch on this batch
        cmd = [
            'vsearch',
            '--usearch_global', str(temp_query),
            '--db', str(self.vsearch_db),
            '--userout', '/dev/stdout',
            '--userfields', 'query+target+id+alnlen+qcov',
            '--id', str(self.min_identity),
            '--query_cov', str(min_coverage),
            '--maxaccepts', str(MAX_TARGET_SEQS),
            '--threads', str(self.threads),
            '--output_no_hits'
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.debug(f"vsearch batch found {len(result.stdout.splitlines())} hits")
        except FileNotFoundError:
            raise RuntimeError("vsearch not found. Please ensure vsearch is installed and in PATH.")
        except subprocess.CalledProcessError as e:
            logging.error(f"vsearch search failed: {e.stderr}")
            raise RuntimeError(f"vsearch search failed: {e.stderr}")

        # Process vsearch results for this batch
        batch_best_hits = {}

        for line in result.stdout.splitlines():
            if not line.strip():
                continue

            parts = line.split('\t')
            if len(parts) != 5:
                continue

            query_id, sseqid, pident, length, coverage = parts
            pident = float(pident)
            length = int(length)
            coverage = float(coverage)
            bitscore = 0.0  # vsearch doesn't compute bitscores

            # Apply minimum length filter (identity and coverage already filtered by vsearch)
            if length >= min_length:
                if sseqid in self.ref_info and self.ref_info[sseqid].genus:
                    ref = self.ref_info[sseqid]
                    # Use genus or species depending on taxonomic_level
                    taxon = ref.species if self.taxonomic_level == "species" else ref.genus
                    # Keep hit if it's better than what we've seen for this read
                    # Since vsearch doesn't have bitscore, use pident as tiebreaker
                    if (query_id not in batch_best_hits or
                            pident > batch_best_hits[query_id][1]):
                        batch_best_hits[query_id] = (taxon, pident, length, coverage, bitscore)
                        if logging.getLogger().isEnabledFor(logging.DEBUG):
                            logging.debug(
                                f"{query_id}: {taxon} "
                                f"(pident={pident:.1f}%, len={length}, "
                                f"cov={coverage:.1f}%)"
                            )

        return batch_best_hits

    def run_search(self, reads_file: Path, min_coverage: float = DEFAULT_MIN_COVERAGE,
                   min_length: int = DEFAULT_MIN_LENGTH) -> ReadSummary:
        """Run sequence search (BLAST or vsearch) on reads file and summarize results in batches.

        Args:
            reads_file: Path to reads file
            min_coverage: Minimum query coverage (0.0-1.0)
            min_length: Minimum alignment length

        Returns:
            ReadSummary with taxonomic counts and per-read details

        Raises:
            ValueError: If database not created or invalid parameters
            RuntimeError: If search operations fail
        """
        if self.search_method == "blast" and not self.blast_db:
            raise ValueError("BLAST database not created. Call create_blast_db() first.")
        if self.search_method == "vsearch" and not self.vsearch_db:
            raise ValueError("vsearch database not created. Call create_vsearch_db() first.")

        # Validate parameters
        if not 0.0 <= min_coverage <= 1.0:
            raise ValueError(f"min_coverage must be between 0.0 and 1.0, got {min_coverage}")
        if min_length < 1:
            raise ValueError(f"min_length must be >= 1, got {min_length}")

        format = self._get_sequence_format(reads_file)

        # Load all reads
        read_records = list(SeqIO.parse(str(reads_file), format))
        total_reads = len(read_records)

        if total_reads == 0:
            return ReadSummary(total_reads=0, taxon_counts=Counter(), unknown_count=0,
                             read_details={}, filtered_low_quality=0, filtered_too_long=0)

        # Apply length filtering
        read_records, num_filtered_length = self._filter_by_length(read_records)
        if num_filtered_length > 0:
            logging.info(f"Filtered {num_filtered_length} reads longer than {self.max_seq_length} bp")

        # Apply quality filtering if enabled
        read_records, num_filtered_quality = self._filter_by_quality(read_records)
        if num_filtered_quality > 0:
            logging.info(f"Filtered {num_filtered_quality} reads below quality threshold {self.min_quality}")

        if len(read_records) == 0:
            logging.warning("All reads filtered by quality/length thresholds")
            return ReadSummary(total_reads=total_reads, taxon_counts=Counter(), unknown_count=0,
                             read_details={}, filtered_low_quality=num_filtered_quality,
                             filtered_too_long=num_filtered_length)

        logging.info(f"Processing {len(read_records)} reads in batches of {self.batch_size}")

        # Process reads in batches
        read_best_hits = {}
        num_batches = (len(read_records) + self.batch_size - 1) // self.batch_size

        # Select batch processing method based on search_method
        batch_method = self.run_blast_batch if self.search_method == "blast" else self.run_vsearch_batch
        desc = f"{self.search_method} batches"

        for batch_num in tqdm(range(num_batches), desc=desc, unit="batch"):
            start_idx = batch_num * self.batch_size
            end_idx = min(start_idx + self.batch_size, len(read_records))
            batch_records = read_records[start_idx:end_idx]

            batch_hits = batch_method(batch_records, min_coverage, min_length)
            read_best_hits.update(batch_hits)

        # Count taxa from best hits
        taxon_counts = Counter()
        for taxon, *_ in read_best_hits.values():
            taxon_counts[taxon] += 1

        # Apply minimum reads filter
        if self.min_reads > 1:
            filtered_taxa = {taxon: count for taxon, count in taxon_counts.items()
                           if count >= self.min_reads}
            num_filtered_taxa = len(taxon_counts) - len(filtered_taxa)
            if num_filtered_taxa > 0:
                logging.info(f"Filtered {num_filtered_taxa} taxa with fewer than {self.min_reads} reads")
            taxon_counts = Counter(filtered_taxa)

        matched_reads = len(read_best_hits)
        unknown_count = len(read_records) - matched_reads + num_filtered_quality + num_filtered_length

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(f"File: {reads_file.name}")
            logging.debug(f"Total reads: {total_reads}")
            logging.debug(f"Matched reads: {matched_reads}")
            logging.debug(f"Unknown reads: {unknown_count}")
            for taxon, count in taxon_counts.most_common():
                logging.debug(f"{taxon}: {count} reads ({count / total_reads * 100:.1f}%)")

        return ReadSummary(total_reads=total_reads,
                           taxon_counts=taxon_counts,
                           unknown_count=unknown_count,
                           read_details=read_best_hits,
                           filtered_low_quality=num_filtered_quality,
                           filtered_too_long=num_filtered_length)

    def run_blast(self, reads_file: Path, min_coverage: float = DEFAULT_MIN_COVERAGE,
                  min_length: int = DEFAULT_MIN_LENGTH) -> ReadSummary:
        """Legacy method: Run BLAST on reads file and summarize results.

        This method is deprecated. Use run_search() instead.
        """
        if self.search_method != "blast":
            raise ValueError("run_blast() can only be used when search_method='blast'")
        return self.run_search(reads_file, min_coverage, min_length)

    def format_summary(self, filename: str, summary: ReadSummary) -> str:
        """Format summary statistics into desired output format.

        Args:
            filename: Name of input file for display (can be str or Path)
            summary: ReadSummary object with results

        Returns:
            Formatted multi-line string with taxonomic breakdown
        """
        # Convert Path to string if needed
        filename_str = str(filename)

        if summary.total_reads == 0:
            return f"{filename_str}\t0\tNo reads found"

        # Sort taxa by count and take top N (or all if max_taxa is None)
        if self.max_taxa is None:
            top_taxa = summary.taxon_counts.most_common()
        else:
            top_taxa = summary.taxon_counts.most_common(self.max_taxa)

        # Calculate percentages and format each taxon on a separate line
        results = []
        header_parts = [filename_str, str(summary.total_reads)]
        filter_notes = []
        if summary.filtered_too_long > 0:
            filter_notes.append(f"{summary.filtered_too_long} too long")
        if summary.filtered_low_quality > 0:
            filter_notes.append(f"{summary.filtered_low_quality} low quality")
        if filter_notes:
            header_parts.append(f"(filtered: {', '.join(filter_notes)})")
        results.append('\t'.join(header_parts))

        other_count = 0
        for taxon, count in top_taxa:
            pct = (count / summary.total_reads) * 100
            if self.max_taxa is None or len(top_taxa) <= self.max_taxa:
                results.append(f"  {pct:.1f}%\t{count}\t{taxon}")
            elif len(results) < self.max_taxa + 1:  # +1 for header line
                results.append(f"  {pct:.1f}%\t{count}\t{taxon}")
            else:
                other_count += count

        # Add Other category if needed (only when max_taxa is set and there are more taxa)
        if self.max_taxa is not None:
            remaining_count = sum(count for taxon, count in summary.taxon_counts.items()
                                  if taxon not in dict(top_taxa))
            other_count += remaining_count
            if other_count > 0:
                other_pct = (other_count / summary.total_reads) * 100
                results.append(f"  {other_pct:.1f}%\t{other_count}\tOther")

        # Add Unknown category if needed
        if summary.unknown_count > 0:
            unknown_pct = (summary.unknown_count / summary.total_reads) * 100
            results.append(f"  {unknown_pct:.1f}%\t{summary.unknown_count}\tUnknown")

        return '\n'.join(results)

    def write_detailed_tsv(self, reads_file: Path, summary: ReadSummary, output_file: Path):
        """Write detailed per-read results to a TSV file.

        Args:
            reads_file: Original reads file path (for filename in output)
            summary: ReadSummary object containing read_details
            output_file: Path to output TSV file
        """
        with open(output_file, 'w') as f:
            # Write header
            f.write("read_id\ttaxon\tpident\talign_length\tcoverage\tbitscore\n")

            # Write matched reads
            for read_id, (taxon, pident, length, coverage, bitscore) in sorted(summary.read_details.items()):
                f.write(f"{read_id}\t{taxon}\t{pident:.2f}\t{length}\t{coverage:.2f}\t{bitscore:.2f}\n")

        logging.info(f"Wrote detailed results to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize read files at taxonomic level against reference sequences"
    )
    parser.add_argument("references", type=Path,
                        help="Reference sequences file (FASTA or FASTQ)")
    parser.add_argument("reads", type=Path, nargs='+',
                        help="Input reads files (FASTA or FASTQ)")
    parser.add_argument("--min-identity", type=float, default=DEFAULT_MIN_IDENTITY,
                        help=f"Minimum similarity threshold for taxon assignment (default: {DEFAULT_MIN_IDENTITY})")
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE,
                        help=f"Minimum query coverage (default: {DEFAULT_MIN_COVERAGE})")
    parser.add_argument("--min-length", type=int, default=DEFAULT_MIN_LENGTH,
                        help=f"Minimum alignment length (default: {DEFAULT_MIN_LENGTH})")
    parser.add_argument("--max-taxa", type=int, default=None,
                        help="Maximum number of taxa to show per file (default: unlimited)")
    parser.add_argument("--threads", type=int, default=-1,
                        help="Number of threads for BLAST searches (default: num cores)")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
                        help=f"Number of queries per BLAST batch (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--level", choices=["genus", "species"], default="species",
                        help="Taxonomic level for reporting (default: species)")
    parser.add_argument("--min-reads", type=int, default=DEFAULT_MIN_READS,
                        help=f"Minimum reads per taxon to report (default: {DEFAULT_MIN_READS})")
    parser.add_argument("--min-quality", type=float, default=None,
                        help="Minimum mean phred quality score (default: no filtering)")
    parser.add_argument("--max-seq-length", type=int, default=DEFAULT_MAX_SEQ_LENGTH,
                        help=f"Maximum sequence length to process (default: {DEFAULT_MAX_SEQ_LENGTH})")
    parser.add_argument("--no-detailed-tsv", action="store_true",
                        help="Don't write detailed per-read results to .tsv files")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress all output except results")
    parser.add_argument("--debug", action="store_true",
                        help="Enable debug logging")
    parser.add_argument("--search-method", choices=["blast", "vsearch"], default="blast",
                        help="Search method to use (default: blast)")

    args = parser.parse_args()

    # Setup logging
    if args.quiet:
        logging.basicConfig(level=logging.ERROR)
    else:
        log_level = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

    if args.threads == -1:
        n_threads = multiprocessing.cpu_count()
    else:
        n_threads = args.threads

    logging.info(f"Will run {args.search_method} with {n_threads} threads")

    # Initialize and run summarizer with context manager for automatic cleanup
    with TaxonSummarizer(
        similarity_threshold=args.min_identity,
        max_taxa=args.max_taxa,
        threads=n_threads,
        batch_size=args.batch_size,
        taxonomic_level=args.level,
        min_reads=args.min_reads,
        min_quality=args.min_quality,
        max_seq_length=args.max_seq_length,
        search_method=args.search_method
    ) as summarizer:
        # Create database from references (BLAST or vsearch)
        if args.search_method == "blast":
            summarizer.create_blast_db(args.references)
        else:
            summarizer.create_vsearch_db(args.references)

        # Process each input file
        for reads_file in args.reads:
            summary = summarizer.run_search(reads_file, args.min_coverage, args.min_length)
            print(summarizer.format_summary(reads_file, summary))

            # Write detailed TSV by default (unless --no-detailed-tsv is specified)
            if not args.no_detailed_tsv:
                # Write to current working directory, not input file's directory
                tsv_filename = reads_file.stem + '.taxonomize.tsv'
                tsv_file = Path.cwd() / tsv_filename
                summarizer.write_detailed_tsv(reads_file, summary, tsv_file)


if __name__ == "__main__":
    main()