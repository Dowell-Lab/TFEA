# Transcription Factor Enrichment Analysis (TFEA)
# Table of Contents
1. <A href="#Help">Help</A>
2. <A href="">Pipeline</A>
1. <A href="#Installation">Installation</A>
2. <A href="#Requirements">Requirements</A>
   - <A href="#Python3">Python3</A>
     - <A href="#HTSeq">HTSeq</A>
     - <A href="#Pybedtools">Pybedtools</A>
     - <A href="#psutil">psutil</A>
   - <A href="#DESeq">DESeq</A>
   - <A href="#Bedtools">Bedtools</A>
   - <A href="#MEMESuite">MEME Suite</A>
     - <A href="#ImageMagick">Image Magick</A>
   - <A href="#FIJIModules">FIJI Modules</A>
3. <A href="#ConfigurationFile">Configuration File</A>
4. <A href="#RunningLocally">Running Locally</A>
5. <A href="#UsingSBATCH">Using SBATCH</A>
6. <A href="#ExampleOutput">Example Output</A>
7. <A href="#ContactInformation">Contact Information</A>
  
<H2 id="Help">Help</H2>

```
usage: TFEA [-h] [--output OUTPUT] [--bed1 [BED1 [BED1 ...]]]
            [--bed2 [BED2 [BED2 ...]]] [--bam1 [BAM1 [BAM1 ...]]]
            [--bam2 [BAM2 [BAM2 ...]]] [--label1 LABEL1] [--label2 LABEL2]
            [--config CONFIG] [--sbatch SBATCH] [--test-install] [--test-full]
            [--combined_file COMBINED_FILE] [--ranked_file RANKED_FILE]
            [--fasta_file FASTA_FILE] [--md] [--mdd]
            [--md_bedfile1 MD_BEDFILE1] [--md_bedfile2 MD_BEDFILE2]
            [--mdd_bedfile1 MDD_BEDFILE1] [--mdd_bedfile2 MDD_BEDFILE2]
            [--md_fasta1 MD_FASTA1] [--md_fasta2 MD_FASTA2]
            [--mdd_fasta1 MDD_FASTA1] [--mdd_fasta2 MDD_FASTA2]
            [--mdd_pval MDD_PVAL] [--mdd_percent MDD_PERCENT]
            [--combine {intersect/merge,merge all,tfit clean,tfit remove small,False}]
            [--rank {deseq,fc,False}] [--scanner {fimo,genome hits}]
            [--enrichment {auc,auc_bgcorrect}] [--debug]
            [--genomefasta GENOMEFASTA] [--fimo_thresh FIMO_THRESH]
            [--fimo_motifs FIMO_MOTIFS] [--fimo_background FIMO_BACKGROUND]
            [--genomehits GENOMEHITS] [--singlemotif SINGLEMOTIF]
            [--permutations PERMUTATIONS] [--largewindow LARGEWINDOW]
            [--smallwindow SMALLWINDOW] [--dpi DPI] [--padjcutoff PADJCUTOFF]
            [--pvalcutoff PVALCUTOFF] [--plotall] [--output_type {txt,html}]
            [--cpus CPUS]

Transcription Factor Enrichment Analysis (TFEA)

optional arguments:
  -h, --help            show this help message and exit

Main Inputs:
  Inputs required for full pipeline

  --output OUTPUT, -o OUTPUT
                        Full path to output directory. If it exists, overwrite
                        its contents.
  --bed1 [BED1 [BED1 ...]]
                        Bed files associated with condition 1
  --bed2 [BED2 [BED2 ...]]
                        Bed files associated with condition 2
  --bam1 [BAM1 [BAM1 ...]]
                        Sorted bam files associated with condition 1. Must be
                        indexed.
  --bam2 [BAM2 [BAM2 ...]]
                        Sorted bam files associated with condition 2. Must be
                        indexed.
  --label1 LABEL1       An informative label for condition 1
  --label2 LABEL2       An informative label for condition 2
  --config CONFIG, -c CONFIG
                        A configuration file that a user may use instead of
                        specifying flags. Command line flags will overwrite
                        options within the config file. See examples in the
                        config_files folder.
  --sbatch SBATCH, -s SBATCH
                        Submits an sbatch (slurm) job. If specified, input an
                        e-mail address.
  --test-install, -ti   Checks whether all requirements are installed and
                        command-line runnable.
  --test-full, -t       Performs unit testing on full TFEA pipeline.

Processed Inputs:
  Input options for pre-processed data

  --combined_file COMBINED_FILE
                        A single bed file combining regions of interest.
  --ranked_file RANKED_FILE
                        A bed file containing each regions rank as the 4th
                        column.
  --fasta_file FASTA_FILE
                        A fasta file containing sequences to be analyzed,
                        ranked by the user.

Secondary Analysis Inputs:
  Input options for performing MD-Score and Differential MD-Score analysis

  --md                  Switch that controls whether to perform MD analysis.
  --mdd                 Switch that controls whether to perform differential
                        MD analysis.
  --md_bedfile1 MD_BEDFILE1
                        A bed file for MD-Score analysis associated with
                        condition 1.
  --md_bedfile2 MD_BEDFILE2
                        A bed file for MD-Score analysis associated with
                        condition 2.
  --mdd_bedfile1 MDD_BEDFILE1
                        A bed file for Differential MD-Score analysis
                        associated with condition 1.
  --mdd_bedfile2 MDD_BEDFILE2
                        A bed file for Differential MD-Score analysis
                        associated with condition 2.
  --md_fasta1 MD_FASTA1
                        A fasta file for MD-Score analysis associated with
                        condition 1.
  --md_fasta2 MD_FASTA2
                        A fasta file for MD-Score analysis associated with
                        condition 2.
  --mdd_fasta1 MDD_FASTA1
                        A fasta file for Differential MD-Score analysis
                        associated with condition 1.
  --mdd_fasta2 MDD_FASTA2
                        A fasta file for Differential MD-Score analysis
                        associated with condition 2.
  --mdd_pval MDD_PVAL   P-value cutoff for retaining differential regions.
                        Default: 0.2
  --mdd_percent MDD_PERCENT
                        Percentage cutoff for retaining differential regions.
                        Default: False

Module Switches:
  Switches for different modules

  --combine {intersect/merge,merge all,tfit clean,tfit remove small,False}
                        Method for combining input bed files
  --rank {deseq,fc,False}
                        Method for ranking combined bed file
  --scanner {fimo,genome hits}
                        Method for scanning fasta files for motifs
  --enrichment {auc,auc_bgcorrect}
                        Method for calculating enrichment
  --debug               Print memory and CPU usage to stderr

Fasta Options:
  Options for performing conversion of bed to fasta

  --genomefasta GENOMEFASTA
                        Genomic fasta file

Scanner Options:
  Options for performing motif scanning

  --fimo_thresh FIMO_THRESH
                        P-value threshold for calling FIMO motif hits.
                        Default: 1e-6
  --fimo_motifs FIMO_MOTIFS
                        Full path to a .meme formatted motif databse file.
                        Some databases included in motif_databases folder.
  --fimo_background FIMO_BACKGROUND
                        Options for choosing mononucleotide background
                        distribution to use with FIMO. {'largewindow',
                        'smallwindow', int, file}
  --genomehits GENOMEHITS
                        A folder containing bed files with pre-calculated
                        motif hits to a genome. For use with 'genome hits'
                        scanner option.
  --singlemotif SINGLEMOTIF
                        Option to run analysis on a subset of motifs within
                        specified motif database or genome hits. Can be a
                        single motif or a comma-separated list of motifs.

Enrichment Options:
  Options for performing enrichment analysis

  --permutations PERMUTATIONS
                        Number of permutations to perfrom for calculating
                        p-value. Default: 1000
  --largewindow LARGEWINDOW
                        The size (bp) of a large window around input regions
                        that captures background. Default: 1500
  --smallwindow SMALLWINDOW
                        The size (bp) of a small window arount input regions
                        that captures signal. Default: 150

Output Options:
  Options for the output.

  --dpi DPI             Resolution of output figures. Default: 100
  --padjcutoff PADJCUTOFF
                        A p-adjusted cutoff value that determines some
                        plotting output.
  --pvalcutoff PVALCUTOFF
                        A p-value cutoff value that determines some plotting
                        output.
  --plotall             Plot graphs for all motifs.Warning: This will make
                        TFEA run much slower andwill result in a large output
                        folder.
  --output_type {txt,html}
                        Specify output type. Selecting html will increase
                        processing time and memory usage. Default: txt

Miscellaneous Options:
  Other options.

  --cpus CPUS           Number of processes to run in parallel. Note:
                        Increasing this value will significantly increase
                        memory footprint. Default: 1
```
 
<H2 id="Pipeline">TFEA Pipeline</H2>
 
![TFEA Pipeline](https://github.com/jdrubin91/TFEA/blob/master/TFEA_Pipeline2.png)
 
<br></br>

<H2 id="Installation">Installation</H2>

To install, this package and all python3 dependencies:

```
git clone https://github.com/jdrubin91/TFEA.git
cd /full/path/to/TFEA/
pip3 install --user .
```

TFEA can then be run from anywhere, try:

```
TFEA --help
```

Alternatively, TFEA can be run without installation using:

```
cd /full/path/to/TFEA/
python3 TFEA/ --help
```
but this will require you to install each required python3 package separately.

It is recommended that you run the multiple test modules within TFEA to make sure everything is properly installed:
```
TFEA --test-install
TFEA --test-full
```

Once you've ran the tests successfully, you should be ready to run the full version of TFEA. Below is the minimum required input to run the full TFEA ppeline. These files are provided within 'test_files' for you to get familiar with TFEA (for all paths, make sure you input the full path to 'test_files' on your machine or ```cd``` into the TFEA parent directory and simply copy paste):

```
#Using only flags
TFEA --output ./test_files/test_rep2 \
--bed1 ./test_files/SRR1105736.tfit_bidirs.chr22.bed ./test_files/SRR1105737.tfit_bidirs.chr22.bed \
--bed2 ./test_files/SRR1105738.tfit_bidirs.chr22.bed ./test_files/SRR1105739.tfit_bidirs.chr22.bed \
--bam1 ./test_files/SRR1105736.sorted.chr22.subsample.bam ./test_files/SRR1105737.sorted.chr22.subsample.bam \
--bam2 ./test_files/SRR1105738.sorted.chr22.subsample.bam ./test_files/SRR1105739.sorted.chr22.subsample.bam \
--label1 condition1 --label2 condition2 \
--genomefasta ./test_files/chr22.fa \
--fimo_motifs ./test_files/test_database.meme

#Using only a config file
TFEA --config ./test_files/test_config.ini

#On FIJI
TFEA --config ./test_files/test_config.ini --sbatch your_email@address.com
```

<br></br>

<H2 id="Requirements">Requirements</H2>

Before running TFEA, make sure you have the following installed on your machine

*Note:* If using the --sbatch module on FIJI, you will need to install DESeq and DESeq2 but you should not need to install anything else
  <H3 id="Python3">Python 3</H3>
  TFEA is written in python3. Instructions for installing python3 can be found here:
  <a href="https://www.python.org/downloads/">python3 Installation</a>

  <H3 id="HTSeq">HTSeq</H3>
  TFEA uses HTSeq to draw meta plots of coverage (from BAM inputs) over inputted BED files.
  
  ```
  pip3 install htseq
  ```
  
  If you are on FIJI compute cluster, htseq is available as a module:
  
  ```
  module load python/3.6.3/htseq
  ```

  <H3 id="DESeq">DESeq</H3>
  TFEA uses DESeq or DESeq2 (depending on replicate number) to rank inputted bed files based on fold change significance. Make sure DESeq and DESeq2 are both installed on your system R, in your terminal:
    
  ```
  R
  > source("https://bioconductor.org/biocLite.R")
  > biocLite("DESeq")
  > biocLite("DESeq2")
  ```
  
  If on FIJI, make sure all gcc modules are unloaded before installing DESeq or DESeq2. This can be accomplished with:
  
  ```
  module unload gcc
  ```
  
  or
  
  ```
  module purge
  ```
  
  <H3 id="Bedtools">Bedtools</H3>
  TFEA uses Bedtools to do several genomic computations. Instructions for installing bedtools can be found here:
  <br></br>
  <a href="http://bedtools.readthedocs.io/en/latest/content/installation.html">Bedtools Installation</a>
  <br></br>
  If you are on FIJI compute cluster, bedtools is available as a module:
  
  ```
  module load bedtools/2.25.0
  ```
  
  <H3 id="MEMESuite">MEME Suite</H3>
  TFEA uses the MEME suite to scan sequences from inputted bed files for motif hits using the background atcg distribution form inputted bed file regions. TFEA also uses the MEME suite to generate motif logos for html display. Instructions for downloading and installing the MEME suite can be found here:
  <br></br>
  <a href="http://meme-suite.org/doc/install.html?man_type=web">MEME Download and Installation</a>
  Note: TFEA uses some 
  <br></br>
  If you are on FIJI compute cluster, the meme suite is available as a module:
  
  ```
  module load meme/5.0.3
  ```
  
  <H3 id="FIJIModules">FIJI Modules</H3>
  Below is a summary of all FIJI modules needed.
  
  ```
  module load python/3.6.3
  module load python/3.6.3/matplotlib/1.5.1
  module load python/3.6.3/scipy/0.17.1
  module load python/3.6.3/numpy/1.14.1
  module load python/3.6.3/htseq/0.9.1
  module load python/3.6.3/pybedtools/0.7.10

  module load bedtools/2.25.0
  module load meme/5.0.3
  ```

<br></br>

<H2 id="ConfigurationFile">Configuration File</H3>
TFEA can be run exclusively through the command line using flags. Alternatively, TFEA can be run using a configuration file (.ini). If both flag inputs and configuration variable inputs are provided, TFEA uses flag inputs preferrentially. Below is an example of a configuration file

  ```bash
  #Example config.ini file for use with TFEA.

  [MODULES]

  #This module combines bed files from BED and merges them using bedtools. If False, it will assume BEDS[0] contains the bed file of interest (must be a sorted bed file). (boolean)
  COMBINE = 'intersect/merge

  #This module performs bedtools multicov which requires bam files and a bed file. It will count reads for each bam file across all regions in the inputted bed file. (boolean)
  COUNT = True

  #This module performs DESeq and then ranks regions based on the p-value obtained from DESeq, if you set this to false, TFEA will look for the DESeq file within your specified output directory. (boolean)
  DESEQ = True

  #This module performs the bulk of the calculation of TFEA and will most likely take the longest. Unless you just want to generate files, this should usually be set to True. (boolean)
  CALCULATE = True

  #Determines whether ES calculations will be run in parallel. This is recommended to speed up the process.       (boolean)
  POOL=True

  #This module draws dots on enrichment scatter plot based on whether they are less than the specified p-value cutoff. (boolean)
  DRAWPVALCUTOFF = False

  #This module allows you to specify whether you want to perform TFEA for all motifs in the specified database or whether you want to do just one motif. If you want to do a single motif, you must specify the exact name of the motif (ex. SINGLEMOTIF = 'HO_SP3_HUMAN.H10MO.B.bed'). (boolean)
  SINGLEMOTIF = False

  [DATA]
  #Full path to where you want the data to be output. TFEA outputs a folder with results. (string)
  OUTPUT = './'

  #A variable that can be used if wanted. This variable can be referenced later on using ${BEDDIR} (optional string)
  BEDDIR = '/full/path/to/BEDS/'

  #A list of full paths to BED files corresponding to all treatments. One or multiple BED files can be used but they MUST be within a list. (list of strings)
  BEDS = [${BEDDIR}+'BEDNAME1.bed',${BEDDIR}+'BEDNAME2.bed']

  #A variable that can be used if wanted. This variable can be referenced later on using ${BAMDIR} (optional string)
  BAMDIR = '/full/path/to/BAMS/'

  #A list of full paths to sorted BAM files corresponding to treatment 1. (list of strings)
  BAM1 = [${BAMDIR}+'CONDITION1_rep1.sorted.bam',${BAMDIR}+'CONDITION1_rep2.sorted.bam']

  #The name of treatment 1. (string)
  LABEL1 = 'Treatment 1'

  #A list of full paths to sorted BAM files corresponding to treatment 2. (list of strings)
  BAM2 = [${BAMDIR}+'CONDITION2_rep1.sorted.bam',${BAMDIR}+'CONDITION2_rep2.sorted.bam']

  #The name of treatment 2. (string)
  LABEL2 = 'Treatment 2'


  [THRESHOLDS]
  #FDR cut off for calling significant hits (float)
  FDRCUTOFF = 0.1

  #P-value cut off for calling significant hits (float)
  PVALCUTOFF = 0.1

  #Corresponds to the furthest motif hits that will be displayed in the enrichment scatter plot. This does not affect results (float)
  LARGEWINDOW = 1500.0

  #Corresponds to the threshold in which a positive or negative hit is called. Changing this parameter will change your results, only change if you have a good reason to do so. (float)
  SMALLWINDOW = 150.0
  ```

<br></br>

<H2 id="RunningLocally">Running Locally</H2>
If you desire to run TFEA on your local machine, make sure you have the required programs installed. In general, Python packages can be installed using pip, others may require additional installation steps.
<br></br>
Once all dependencies are installed, TFEA can be run using:

  ```bash
  python3 TFEA/ --config CONFIG.ini
  ```

<br></br>

<H2 id="UsingSBATCH">Using SBATCH</H2>
Submitting jobs through the slurm scheduler is supported. To use this module:

  ```bash
  python3 TFEA/ --config CONFIG.ini --sbatch email@address.com
  ```


Node configuration can be changed within scripts/run_main.sbatch. See here the sbatch code used:

  ```qsub
  #!/bin/bash

  ###Name the job
  #SBATCH --job-name=TFEA

  ###Specify the queue
  #SBATCH -p short

  ###Specify WallTime
  #SBATCH --time=24:00:00

  ### Specify the number of nodes/cores
  #SBATCH --ntasks=10

  ### Allocate the amount of memory needed
  #SBATCH --mem=20gb

  ### Set error and output locations. These will be automatically updated to the output directory.
  #SBATCH --error /scratch/Users/user/e_and_o/%x.err
  #SBATCH --output /scratch/Users/user/e_and_o/%x.out

  ### Set your email address. This is changed automatically
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=jonathan.rubin@colorado.edu

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

  python3 ${cmd}
  ```
**NOTE:** For TFEA to properly run a job, the python call within the sbatch script:
>python3 ${cmd}

**MUST NOT BE CHANGED**

<br></br>

<H2 id="ExampleOutput">Example Output</H2>
This part is still under construction.

<br></br>

<H2 id="ContactInformation">Contact Information</H2>
<a href="mailto:jonathan.rubin@colorado.edu">Jonathan.Rubin@colorado.edu</a>

