#!/bin/bash

###Name the job
#SBATCH --job-name=TFEA

###Specify the queue
#SBATCH -p long

###Specify WallTime
#SBATCH --time=10:00:00

### Specify the number of nodes/cores
#SBATCH --nodes=1
#SBATCH --ntasks=1

### Allocate the amount of memory needed
#SBATCH --mem=64gb

### Setting to mail when the job is complete
#SBATCH --error /scratch/Users/joru1876/e_and_o/%x.err
#SBATCH --output /scratch/Users/joru1876/e_and_o/%x.out

### Set your email address
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joru1876@colorado.edu

### Load required modules
module purge
module load python/2.7.14
module load bedtools/2.25.0
module load python/2.7.14/matplotlib/1.5.1
module load python/2.7.14/scipy/0.17.1

### now call your program

src=/scratch/Users/joru1876/TFEA/src/

python $src ${output} ${filedir} ${figuredir} ${e_and_o}
