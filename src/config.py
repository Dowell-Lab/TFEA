#User defined input. Spceify here the parameters used to run TFEA including full paths to all the necessary input files
#======================================================================================================================#
#Which parts of TFEA would you like to run? These are switches to turn on/off different modules in TFEA
COMBINE = False
COUNT = False
RANK = True
DISTANCE = True
CALCULATE = True

#Input a list of bed files with regions of interest to be analyzed. Ideally, these are meant to be Tfit output files corresponding to each bam file submitted. 
#These files will be concatenated and merged (bedtools) to produce detected regions in all samples. If you only have one bed file with regions of interest
#submit it as a single item in the BEDS list and set COMBINE to False.
BEDS=['/full/path/to/bed/file1.bed','/full/path/to/bed/file2.bed','/full/path/to/bed/file3.bed']

#Input bam files as a list containing transcription data. Must specify at least two bam files. 
#If multiple replicates, specify each as a full path in the appropriate list.
BAMDIR = '/Users/joru1876/scratch_backup/TFEA_grant_figure/'
BAM1 = ['/full/path/to/bam1/rep1.bam','/full/path/to/bam1/rep2.bam','...']
BAM2 = ['/full/path/to/bam2/rep1.bam','/full/path/to/bam2/rep2.bam','...']

#Specify conditions for bam files
LABEL1 = 'Condition1'
LABEL2 = 'Condition2'

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
DESEQFILE = '/Users/joru1876/scratch_backup/TFEA_grant_figure/Tfit_SRR3739_bidir_predictions.sorted.merge.count.bed.id.bed.DMSONutlinnascent.res.txt'
