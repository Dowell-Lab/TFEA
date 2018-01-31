# TFEA
Transcription Factor Enrichment Analysis

This repo has the following requirements:
  1. module load python/2.7.14
  2. module load bedtools/2.25.0
  3. module load python/2.7.14/matplotlib/1.5.1
  4. module load python/2.7.14/scipy/0.17.1

The following inputs are required and can be modified in src/config.py
  1. BEDS: A list of full paths to bed files (must be a list) containing regions of interest with which to perform TFEA. Ideally, these are Tfit regions. One bed file can be specified but it still must be within a list.
  2. BAM1: A list of bam files (must be a list) containing all replicates for a condition
  3. BAM2: A list of bam files (must be a list) containing all replicates for another condition
  4. LABEL1: Name of condition one corresponding to BAM1 files
  5. LABEL2: Name of condition two corresponding to BAM2 files
  6. OUTPUT: Full path to output file
  7. CONDITION: Boolean, dictates whether you would like TFEA to get bam files located in pubgro from a conditions table. Set to false if you know the full paths to all your BAMS and BEDS
  8. SINGLEMOTIF: Use if you want to only perform TFEA on a single motif within one of the databases in MEME. If used, must correspond to a motif in the database verbatim
  9. MOTIF_HITS: A directory with motif hits across the genome for all motifs in HOCOMCO (currently this directory has v10 motifs)

For users running this on a compute cluster (i.e. CU Boulder FIJI), it is recommended to run src/sbatch.py which will automatically create necessary output folders and submit the appropriate sbatch script to the slurm scheduler.
