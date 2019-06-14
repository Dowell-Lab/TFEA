#!/bin/bash

###Name the job
#SBATCH --job-name=TFEA_simulate

###Specify the queue
#SBATCH -p long

###Specify WallTime
#SBATCH --time=240:00:00

### Specify the number of nodes/cores
#SBATCH --ntasks=10

### Allocate the amount of memory needed
#SBATCH --mem=50gb

### Setting to mail when the job is complete
#SBATCH --error /scratch/Users/joru1876/e_and_o/%x.err
#SBATCH --output /scratch/Users/joru1876/e_and_o/%x.out

### Set your email address
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joru1876@colorado.edu

### Load required modules
module purge

module load python/3.6.3
module load python/3.6.3/matplotlib/1.5.1
module load python/3.6.3/scipy/0.17.1
module load python/3.6.3/numpy/1.14.1
module load python/3.6.3/htseq/0.9.1
module load python/3.6.3/pybedtools/0.7.10

module load bedtools/2.25.0
module load meme/5.0.3

### now call your program

### python ${src} --config ${config} --sbatch SUBMITTED
python3 ${cmd}
