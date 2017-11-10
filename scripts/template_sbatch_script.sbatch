#!/bin/bash

###Name the job
#SBATCH --job-name=TFEA
### Specify the number of nodes/cores
#SBATCH --nodes=1
#SBATCH --ntasks=1

### Allocate the amount of memory needed
#SBATCH --mem=10gb

### Set your expected walltime
#SBATCH --time=100:00:00

### Setting to mail when the job is complete
#SBATCH --error /Users/joru1876/qsub_errors/
#SBATCH --output /Users/joru1876/qsub_stdo/  

### Set your email address
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joru1876@colorado.edu

### Load required modules
module load bedtools/2.25.0
module load meme/4.12.0

### now call your program

src=/Users/joru1876/scratch_backup/TFEA/src/

python $src

