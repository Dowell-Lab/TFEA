# TFEA
Transcription Factor Enrichment Analysis

This repo has the following requirements:
  1. MEME - version 4.10.2
  2. Bedtools - version 2.22.0

The following inputs are required and can be modified in src/config.py
  1. BED: Full path to a bed file containing regions of interest with which to perform TFEA
  2. BAM1: A list of bam files (must be a list) containing all replicates for a condition
  3. BAM2: A list of bam files (must be a list) containing all replicates for another condition
  4. LABEL1: Name of condition one corresponding to BAM1 files
  5. LABEL2: Name of condition two corresponding to BAM2 files
  6. SINGLEMOTIF: Use if you want to only perform TFEA on a single motif within one of the databases in MEME. If used, must correspond to a motif in the database verbatim
  7. DATABASE: Specify which database to get motif PSSMs from. Only databases in MEME are currently available
  8. GENOME: Full path to a genome fasta file from which to convert your regions of interest into sequences to search for motifs
  9. MEMEDB: Path to where motif databases are located within your MEME install directory
