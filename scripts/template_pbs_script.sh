###Name the job
#PBS -N TFEA
### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=1

### Allocate the amount of memory needed
#PBS -l mem=10gb

### Set your expected walltime
#PBS -l walltime=100:00:00

### Setting to mail when the job is complete
#PBS -e /Users/joru1876/qsub_errors/
#PBS -o /Users/joru1876/qsub_stdo/  

### Set your email address
#PBS -m ae
#PBS -M joru1876@colorado.edu



### Choose your shell 
#PBS -S /bin/bash
### Pass enviroment variables to the job

module load bedtools2_2.22.0
module load meme_4.10.1_4

### now call your program

src=/Users/joru1876/scratch_backup/TFEA/src/

python $src

