import read_conditions

#User defined input. Spceify here the parameters used to run TFEA including full paths to all the necessary input files
#======================================================================================================================#
#Will you be using the batch analysis module? If so, BEDS, BAM variables, and LABEL variables are specified in sbatch script.
BATCH = False

#Which parts of TFEA would you like to run? These are switches to turn on/off different modules in TFEA
COMBINE = True
COUNT = True
RANK = True
DISTANCE = True
CALCULATE = True

#Specify output directory
OUTPUT = False

#Define the FDR cutoff to be used for calling significant hits
FDRCUTOFF = pow(10,-6)
# FDRCUTOFF = 0.05

CONDITION = False

if CONDITION:
    CONDITIONS='/scratch/Users/joru1876/TFEA_files/conditions_short_20161103_tentative.txt_20161107-165140.csv'
    SPECIFICCELLTYPE = 'HCT116'
    LABEL1='DMSO_1hr'
    LABEL2='Nutlin_1hr'
    BAMDIR='/scratch/Users/joru1876/TFEA_files/bams/'
    BEDDIR='/scratch/Users/joru1876/TFEA_files/tfit_beds/'
    KEYWORD='Allen2014'
    FILEDIR='/scratch/Users/joru1876/TFEA_files/Allen2014-30/'
    BAM1,BAM2,BEDS = read_conditions.run(CONDITIONS,KEYWORD,SPECIFICCELLTYPE,LABEL1,LABEL2,BAMDIR,BEDDIR)
else:
    #Input a list of bed files with regions of interest to be analyzed. Ideally, these are meant to be Tfit output files corresponding to each bam file submitted. 
    #These files will be concatenated and merged (bedtools) to produce detected regions in all samples. If you only have one bed file with regions of interest
    #submit it as a single item in the BEDS list and set COMBINE to False.
    # BEDDIR1 = ''
    BEDDIR2 = '/scratch/Users/joru1876/Taatjes/171026_NB501447_0180_fastq_IRISREP2/Demux/Taatjes-374/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit/'
    BEDS = [BEDDIR2+'JDR_Tfit-5_bidir_predictions.bed',BEDDIR2+'JDR_Tfit-7_bidir_predictions.bed']

    #Input bam files as a list containing transcription data. Must specify at least two bam files. 
    #If multiple replicates, specify each as a full path in the appropriate list.
    BAMDIR1 = '/scratch/Users/joru1876/Taatjes/161220_K00262_0062_BHH7CHBBXX_IRISREP1/trimmed/flipped/bowtie2_first_run/sortedbam/'
    BAMDIR2 = '/scratch/Users/joru1876/Taatjes/171026_NB501447_0180_fastq_IRISREP2/Demux/Taatjes-374/trimmed/flipped/bowtie2/sortedbam/'
    BAM1 = [BAMDIR2+'0_2_S1_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam',BAMDIR1+'0-1_S11_L006_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam']
    BAM2 = [BAMDIR2+'30_2_S3_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam',BAMDIR1+'30-1_S15_L007_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam']

    #Specify conditions for bam files
    LABEL1 = '0 IFN'
    LABEL2 = '30 IFN'

    FILEDIR = '/scratch/Users/joru1876/TFEA_files/IRIS/'

#Specify cell type
CELLTYPE = 'HCT116'

#Specify whether you want to run a single motif or a database. By default, TFEA runs on the latest version of HOCMOCO obtained through MEME.
#Default:False. Change to a motif name if you want to run TFEA on a single motif (make sure your single motif is in the specified database).
# SINGLEMOTIF='HO_P53_HUMAN.H10MO.B.bed'
SINGLEMOTIF=False
DATABASE='HOCOMOCOv11_full_HUMAN_mono_meme_format.meme'

#Specify full path to genome fasta file
GENOME = '/Users/joru1876/scratch_backup/hg19_reference_files/hg19_all.fa'

#Specify full path to MEME directory with motif databases. Make sure you have updated the motif databases in MEME 
#using the following command (where MEMEDB is the full path to where MEME databases are located):
#update-sequence-db MEMEDB
MEMEDB = '/scratch/Users/joru1876/scratch_backup/TFEA/motif_databases/'

MOTIF_HITS = '/scratch/Shares/dowell/md_score_paper/PSSM_hits_genome_wide/pval-6/'

#Optional: If choosing deseqfile() option for ranking, provide deseq res.txt output file 
# DESEQFILE = '/Users/joru1876/scratch_backup/TFEA_grant_figure/Tfit_SRR3739_bidir_predictions.sorted.merge.count.bed.id.bed.DMSONutlinnascent.res.txt'
DESEQFILE = '/scratch/Users/joru1876/TFEA_files/Allen2014/DESeq.res.txt'
