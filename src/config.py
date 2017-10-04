#User defined input. Spceify here the parameters used to run TFEA including full paths to all the necessary input files
#======================================================================================================================#
#Input bed file containing regions of interest to perform TFEA. This file is required to run TFEA. Ideally a list of eRNAs called by Tfit.
BED='full/path/to/bed/file.bed'

#Input bam files as a list containing transcription data. Must specify at least two bam files. If multiple replicates, specify each as a full path in the appropriate list.
BAM1 = ['full/path/to/bam1/rep1.bam','full/path/to/bam1/rep2.bam','...']
BAM2 = ['full/path/to/bam2/rep1.bam','full/path/to/bam2/rep2.bam','...']

#Specify conditions for bam files
LABEL1 = 'Condition1'
LABEL2 = 'Condition2'

#Specify whether you want to run a single motif or a database. By default, TFEA runs on the latest version of HOCMOCO obtained through MEME.
#Default:False. Change to a motif name if you want to run TFEA on a single motif (make sure your single motif is in the specified database).
SINGLEMOTIF=False
DATABASE='hocomoco'

#Specify full path to genome fasta file
GENOME = 'full/path/to/genome.fa'