import read_conditions

#User defined input. Spceify here the parameters used to run TFEA including full paths to all the necessary input files
#======================================================================================================================#

#Which parts of TFEA would you like to run? These are switches to turn on/off different modules in TFEA

#This module combines bed files from BED and merges them using bedtools. If False, it will assume BEDS[0] contains the bed file of interest (must be a sorted bed file)
COMBINE = True

#This module performs bedtools multicov which requires bam files and a bed file. It will count reads for each bam file across all regions in the inputted bed file
COUNT = True

#This module performs DESeq and then ranks regions based on the p-value obtained from DESeq, if you set this to false, TFEA will look for the DESeq file within your specified output directory
DESEQ = True

#This module performs the bulk of the calculation of TFEA and will most likely take the longest. Unless you just want to generate files, this should usually be set to True
CALCULATE = True

#Define the FDR cutoff to be used for calling significant hits
# FDRCUTOFF = pow(10,-6)
FDRCUTOFF = 0.05
PVALCUTOFF = 0.01
DRAWPVALCUTOFF = False

CONDITION = False

if CONDITION:
    CONDITIONS='/scratch/Users/joru1876/TFEA_files/conditions_short_20161103_tentative.txt_20161107-165140.csv'
    SPECIFICCELLTYPE = 'MCF7'
    LABEL1='vehicle'
    LABEL2='E2_40min'
    BAMDIR='/scratch/Users/joru1876/TFEA_files/bams/'
    BEDDIR='/scratch/Users/joru1876/TFEA_files/tfit_beds/'
    KEYWORD='Hah2013'
    OUTPUT='/scratch/Users/joru1876/TFEA_files/Hah2013/'
    BAM1,BAM2,BEDS = read_conditions.run(CONDITIONS,KEYWORD,SPECIFICCELLTYPE,LABEL1,LABEL2,BAMDIR,BEDDIR)
else:
    # #Input a list of bed files with regions of interest to be analyzed. Ideally, these are meant to be Tfit output files corresponding to each bam file submitted. 
    # #These files will be concatenated and merged (bedtools) to produce detected regions in all samples. If you only have one bed file with regions of interest
    # #submit it as a single item in the BEDS list and set COMBINE to False.
    # BEDDIR2 = '/scratch/Users/joru1876/Taatjes/171026_NB501447_0180_fastq_IRISREP2/Demux/Taatjes-374/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit/'
    # BEDS = [BEDDIR2+'JDR_Tfit-5_bidir_predictions.bed',BEDDIR2+'JDR_Tfit-7_bidir_predictions.bed']

    # #Input bam files as a list containing transcription data. Must specify at least two bam files. 
    # #If multiple replicates, specify each as a full path in the appropriate list.
    # BAMDIR1 = '/scratch/Users/joru1876/Taatjes/161220_K00262_0062_BHH7CHBBXX_IRISREP1/trimmed/flipped/bowtie2_first_run/sortedbam/'
    # BAMDIR2 = '/scratch/Users/joru1876/Taatjes/171026_NB501447_0180_fastq_IRISREP2/Demux/Taatjes-374/trimmed/flipped/bowtie2/sortedbam/'
    # BAM1 = [BAMDIR2+'0_2_S1_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam',BAMDIR1+'0-1_S11_L006_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam']
    # BAM2 = [BAMDIR2+'30_2_S3_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam',BAMDIR1+'30-1_S15_L007_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam']

    # #Specify conditions for bam files
    # LABEL1 = '0_IFN'
    # LABEL2 = '30_IFN'

    # OUTPUT = '/scratch/Users/joru1876/TFEA_files/IRIS/'

    #Input a list of bed files with regions of interest to be analyzed. Ideally, these are meant to be Tfit output files corresponding to each bam file submitted. 
    #These files will be concatenated and merged (bedtools) to produce detected regions in all samples. If you only have one bed file with regions of interest
    #submit it as a single item in the BEDS list and set COMBINE to False.
    BEDDIR1 = '/scratch/Users/joru1876/Taatjes/170825_NB501447_0152_fastq_SERCAREP2_30REP1_RESEQUENCING/Demux/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit/'
    BEDDIR2 = '/scratch/Users/joru1876/Taatjes/170207_K00262_0069_AHHMHVBBXX_SERCAREP1/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit_run2/'
    BEDS = [BEDDIR1+'foot_print_testing-9_bidir_predictions.bed',BEDDIR1+'foot_print_testing-11_bidir_predictions.bed',BEDDIR2+'foot_print_testing-7_bidir_predictions.bed',BEDDIR2+'foot_print_testing-11_bidir_predictions.bed']

    #Input bam files as a list containing transcription data. Must specify at least two bam files. 
    #If multiple replicates, specify each as a full path in the appropriate list.
    BAMDIR1 = '/scratch/Users/joru1876/Taatjes/170207_K00262_0069_AHHMHVBBXX_SERCAREP1/cat/trimmed/flipped/bowtie2/sortedbam/'
    BAMDIR2 = '/scratch/Users/joru1876/Taatjes/170825_NB501447_0152_fastq_SERCAREP2_30REP1_RESEQUENCING/Demux/cat/trimmed/flipped/bowtie2/sortedbam/'
    BAM1 = [BAMDIR1+'J12_trimmed.flip.fastq.bowtie2.sorted.bam',BAMDIR2+'J1DO1_AGTCAA_S1_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam']
    BAM2 = [BAMDIR1+'J52_trimmed.flip.fastq.bowtie2.sorted.bam',BAMDIR2+'J5D451_GTCCGC_S3_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam']

    #Specify conditions for bam files (no spaces allwed)
    LABEL1 = '0_Serum'
    LABEL2 = '45_Serum'

    OUTPUT = '/scratch/Users/joru1876/TFEA_files/Rubin/'

#Specify whether you want to run a single motif or a database. By default, TFEA runs on the latest version of HOCMOCO obtained through MEME.
#Default:False. Change to a motif name if you want to run TFEA on a single motif (make sure your single motif is in the specified database).
# SINGLEMOTIF='HO_P53_HUMAN.H10MO.B.bed'
# SINGLEMOTIF='HO_PROX1_HUMAN.H10MO.D.bed'
# SINGLEMOTIF='HO_ZBED1_HUMAN.H10MO.D.bed'
SINGLEMOTIF=False

#This is a folder that contains PSSM hits across the genome. This folder is needed for running DE-Seq and must be downloaded separately
MOTIF_HITS = '/scratch/Shares/dowell/md_score_paper/PSSM_hits_genome_wide/pval-6/'
