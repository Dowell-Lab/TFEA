#!/bin/bash
#SBATCH --job-name=TFEA_test
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=hoto7260@colorado.edu
#SBATCH --ntasks=5
#SBATCH -N 1
#SBATCH -p short 
#SBATCH --mem=10gb
#SBATCH --time=00:30:00 
#SBATCH --output=/Users/hoto7260/src/TFEA/test/out/e_and_o/TFEA_test_%j.out
#SBATCH --error=/Users/hoto7260/src/TFEA/test/out/e_and_o/TFEA_test_%j.out

########################################
########## EDIT THIS ##################
SRC_DIR=~/src/TFEA

####################################################
####  ENSURE THESE POINT TO YOUR ENVIRONMENTS  #####
####################################################
### Clear modules and load conda environment
module purge

module load python/3.7.4
module load samtools/1.3.1
module load bedtools/2.25.0
module load meme/5.0.3
module load samtools/1.3.1
module load gcc/7.1.0
module load R/4.4.0

echo PATH: $PATH
source ${SRC_DIR}/vir_TFEA/bin/activate
echo PATH: $PATH


#################################################
############  SET UP VARIABLES   #################
#################################################
#dir=/scratch/Users/dara6367/PRO-seq_interspecies-nutlin/Bidirectional-Flow
outdir=${SRC_DIR}/test/out
fimo_motifs=${SRC_DIR}/motif_files/H12P53_meme_format.meme
fimo_pval=${SRC_DIR}/assets/PVAL_CUTOFF_ESTIMATES_H12.txt 
fimo_back=${SRC_DIR}/assets/enhancer_background_flat
genomefasta=/scratch/Shares/dowell/genomes/hg38/hg38.fa
input_dir=${SRC_DIR}/TFEA/test/test_files/

###Display job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors, across $SLURM_NNODES nodes, with $SLURM_JOB_CPUS_PER_NODE cpus per node



# make the out directory
mkdir ${outdir}

#cd ${SRC_DIR}
#python setup.py build
#python setup.py install

# #######################################
# ###########  RUN TFEA  ################
# #######################################
echo "Running TFEA with just a consensus file"
#python3 ${SRC_DIR}/vir_TFEA/lib/python3.7/site-packages/tfea-1.1.4-py3.7.egg/TFEA \
#    --output ${outdir}/via_combined_file/ \
#    --combined_file ${input_dir}/test_combined_file.bed \
#    --bam1 ${input_dir}/SRR1105736.sorted.chr22.subsample.bam ${input_dir}/SRR1105737.sorted.chr22.subsample.bam \
#    --bam2 ${input_dir}/SRR1105738.sorted.chr22.subsample.bam ${input_dir}/SRR1105739.sorted.chr22.subsample.bam \
#    --label1 condition1 --label2 condition2 \
#    --genomefasta ${input_dir}/chr22.fa \
#    --fimo_motifs ${input_dir}/test_database.meme \
#    --cpus 5 \
#    --mem 10gb \


echo "Running TFEA with counts (Suggested for Nascent)"
python3 ${SRC_DIR}/vir_TFEA/lib/python3.7/site-packages/tfea-1.1.4-py3.7.egg/TFEA  \
        --output ${outdir}/via_counts2/\
        --count_file ${input_dir}/New_TFEA/fixed_MU_nongenetss_SJSAp53test_8.20.25_counts.txt \
        --label1 VEH \
        --label2 Nutlin \
        --cpus 5 \
        --mem 10gb \
        --treatment VEH,VEH,Nutlin,Nutlin \
        --fimo_motifs ${fimo_motifs} \
        --genomefasta ${genomefasta} \
        --fimo_background ${fimo_back} \
        --fimo_thresh ${SRC_DIR}/assets/PVAL_CUTOFF_ESTIMATES_H12.CORE.txt \
        --output_type html


echo "Running TFEA with a ranked file"
python3 ${SRC_DIR}/vir_TFEA/lib/python3.7/site-packages/tfea-1.1.4-py3.7.egg/TFEA  \
        --output ${outdir}/via_ranked_file2/ \
        --ranked_file /Users/hoto7260/src/TFEA/test_set/ranked_file.bed \
        --label1 VEH \
        --label2 Nutlin \
        --genomefasta ${genomefasta} \
        --fimo_motifs ${fimo_motifs}\
        --cpus 5 \
        --mem 10gb \
        --fimo_background ${fimo_back} \
        --fimo_thresh ${SRC_DIR}/assets/PVAL_CUTOFF_ESTIMATES_H12.CORE.txt \
        --output_type html