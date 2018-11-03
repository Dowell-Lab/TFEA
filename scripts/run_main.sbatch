#!/bin/bash

###Name the job
#SBATCH --job-name=TFEA

###Specify the queue
#SBATCH -p short

###Specify WallTime
#SBATCH --time=4:00:00

### Specify the number of nodes/cores
#SBATCH --nodes=1
#SBATCH --ntasks=10

### Allocate the amount of memory needed
#SBATCH --mem=50gb

### Setting to mail when the job is complete
#SBATCH --error /scratch/Users/rusi2317/e_and_o/%x.err
#SBATCH --output /scratch/Users/rusi2317/e_and_o/%x.out

### Set your email address
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rutendo.sigauke@ucdenver.edu

### Load required modules
module purge
module load python/2.7.14
module load bedtools/2.25.0
module load python/2.7.14/matplotlib/1.5.1
module load python/2.7.14/scipy/0.17.1
module load python/2.7.14/htseq
module load python/2.7.14/seaborn
module load samtools/1.3.1
module load meme/4.12.0
module load python/2.7.14/pandas/0.18.1

### now call your program

### python ${src} --config ${config} --sbatch SUBMITTED
python ${cmd}
