#!/bin/bash
#SBATCH --job-name=TFEA
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=hoto7260@colorado.edu
#SBATCH --ntasks=50
#SBATCH -N 1
#SBATCH -p short 
#SBATCH --mem=50gb
#SBATCH --time=24:00:00 
#SBATCH --output=/scratch/Users/hoto7260/Resp_Env/Sasse2019_3p/e_and_o/TFEA_Bench_%j.out
#SBATCH --error=/scratch/Users/hoto7260/Resp_Env/Sasse2019_3p/e_and_o/TFEA_Bench_%j.out


#################################################
############  SET UP VARIABLES  #################
#################################################
outdir=/scratch/Users/hoto7260/Resp_Env/Sasse2019_3p/TFEA_Bench/out/Norm_12.19_HOCO12/
fimo_motifs=/scratch/Users/hoto7260/assets/H12CORE_meme_format.meme
thresh_file=/scratch/Users/hoto7260/assets/PVAL_CUTOFF_ESTIMATES_H12.CORE.txt
genomefasta=/scratch/Shares/dowell/genomes/hg38/hg38.fa
input_files=(/scratch/Users/hoto7260/Resp_Env/Sasse2019_3p/TFEA_Bench/data/Norm_12.19/testing_iter_LE/*)


###Display job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors, across $SLURM_NNODES nodes, with $SLURM_JOB_CPUS_PER_NODE cpus per node

#################################################
###########  SET UP ENVIRONMENT  ################
#################################################
### Clear modules and load conda environment
module purge

module load python/3.7.4
module load samtools/1.3.1
module load bedtools/2.25.0
module load meme/5.0.3
module load samtools/1.3.1
module load gcc/7.1.0
module load R/4.4.0

# Activate environment for the additional python related items
source /Users/hoto7260/src/TFEA/vir_TFEA/bin/activate

echo "STARTING"


# #######################################
# ###########  RUN TFEA  ################
# #######################################
for ranked_file in "${input_files[@]}"; do
    # Get the name for the ranked file (1. remove prefix ending in /, 2. remove suffix starting with .bed)
    tmp=${ranked_file##*/}
    name=${tmp%.bed*}
    echo DIRECTORY ${outdir}/${name}
    printf "\n################################${name}######################################\n"
    # make the out directory
    #TFEA \
    python3 /Users/hoto7260/src/TFEA/vir_TFEA/lib/python3.7/site-packages/tfea-1.1.4-py3.7.egg/TFEA  \
        --output ${outdir}/${name} \
        --ranked_file ${ranked_file} \
        --label1 ETH \
        --label2 TNF_DEX \
        --genomefasta ${genomefasta} \
        --fimo_motifs ${fimo_motifs} \
        --fimo_thresh ${thresh_file} \
        --fimo_background /scratch/Users/hoto7260/assets/enhancer_background_flat \
        --cpus 50 \
        --mem 50gb \
        --output_type html \
        --debug
    wait
    echo "Previous job completed."
    date
    printf "\n######################################################################\n"
done


date
