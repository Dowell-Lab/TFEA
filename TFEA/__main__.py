#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This is the main script that runs when the TFEA package is run.
'''
#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'
__version__ = '4.0'

#Imports
#==============================================================================
import sys
import subprocess
from pathlib import Path

import process_inputs

#ARGUMENT PARSING
#==============================================================================
'''We begin by parsing user arguments. TFEA can be run in two ways and these
    are not mutually exclusive. TFEA has traditional command line flags that
    a user may specify. Additionally, a user may provide a configuration file
    (.ini) with all necessary inputs. Finally, a user may provide both a 
    configuration file and command line flags. In this case, the command line
    flags will overwrite any redundant options in the configuration file.
'''
#TFEA source directory
srcdirectory = Path(__file__).absolute().parent

#argparse to add arguments to this python package
parser = process_inputs.read_arguments()

#Display help message when no args are passed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

#If user provided arguments, then parse them
sbatch = parser.parse_args().SBATCH
test = parser.parse_args().TEST

#TEST module
#==============================================================================
'''If test flag specified, run unittests and exit.
'''
if test:
    subprocess.call(["python3", "-m", "unittest", "-v", "-f", 
                    srcdirectory / 'test' / 'test_basic.py'])
    sys.exit()

#VERIFICATION OF USER INPUTS
#==============================================================================
'''This section of the code reads config file and user specified flags, makes
    sure these are complete and not conflicting and writes them to config.py
    within TFEA for global use across modules
'''
#Verify inputs to make sure user has all necessary variables to properly run TFEA
process_inputs.verify_arguments(parser=parser)

#CREATING DIRECTORIES
#==============================================================================
'''TFEA creates the specified output directory if it doesn't exist. Within the
    output directory, 3 directories are created: 'temp_files', 'e_and_o', and
    'plots'. These contain temporary files, stderr and stdout files, and
    figures generated by TFEA./
'''
#If user specifies the --sbatch flag, then we first create the output 
#directories then run the sbatch script with the 'SUBMITTED' command submitted 
#to the --sbatch flag so we know not to remake output directories. If --sbatch 
#flag not specified, simply make output directories and continue.
process_inputs.create_directories(srcdirectory=srcdirectory)



#==============================================================================
#MAIN SCRIPT
#==============================================================================

#SECONDARY IMPORTS
#==============================================================================
import config

#Print starting statements
#==============================================================================
print("TFEA start: ", file=sys.stderr)

#Print multiprocessing information to stderr
if config.vars.DEBUG:
    mp.log_to_stderr()
    multiprocess.current_mem_usage()

#COMBINE module
#==============================================================================
'''This module is a pre-processing step where a user may specify how to handle
    multiple bed file inputs. The goal is to arrive at a single bed file to
    input into subsequent modules.
'''
if config.vars.COMBINE != False:
    import combine
    combine.main()

#RANK module
#==============================================================================
'''This module decides how to rank regions within the bed files. If genome
    hits specified then the ranked output will only contain the center of each
    region (since we will perform bedtools closest later)
'''
if config.vars.RANK != False:
    import rank
    rank.main()

#SCANNER module
#==============================================================================
'''This module returns motif distances to regions of interest. This is
    accomplished either by scanning regions on the fly using fimo or homer, or 
    by running bedtools closest on region centers compared to a database of
    motif hits across the genome.
'''
import scanner
scanner.main()
    
#ENRICHMENT module
#==============================================================================
'''Where the bulk of TFEA analysis occurs. Some components of plotting module 
    are contained within this enrichment module
'''
import enrichment
enrichment.main()
    
#OUTPUT module
#==============================================================================
'''A module to write output to either a txt or html file
'''
import output
output.main()

print("TFEA done.", file=sys.stderr)

#REMOVE TEMP FILES
#==============================================================================
# if not config.TEMP:
#     print("Removing temporary bed files...")
#     files_to_keep = ['count_file.bed', 'DESeq.R', 'DESeq.res.txt', 
#                     'DESeq.Rout', 'ranked_file.bed', 'markov_background.txt']
#     for file1 in os.listdir(tempdir):
#         if file1 not in files_to_keep:
#              os.remove(os.path.join(tempdir, file1))

#     print("\t", "done")

#==============================================================================
#END OF MAIN SCRIPT
#==============================================================================