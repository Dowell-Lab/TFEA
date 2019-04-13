#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This file contains a list of independent functions that do not call other 
    functions
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
import shutil
import argparse
import subprocess
from pathlib import Path

import exceptions

#Classes
#==============================================================================
class Bunch(object):
    '''This is a 'Bunch' object. It converts a dict to a namespace representation.
        Created by combining code written by Alex Martelli 
        (https://stackoverflow.com/questions/2597278/python-load-variables-in-a-dict-into-namespace)
        and the built-in python types.SimpleNamespace class.
    '''
    def __init__(self, adict):
        self.__dict__.update(adict)
    def __repr__(self):
        keys = sorted(self.__dict__)
        items = ("{}={!r}".format(k, self.__dict__[k]) for k in keys)
        return "{}({})".format(type(self).__name__, ", ".join(items))

#Functions
#==============================================================================
def read_arguments():
    '''This function consolidates all argument parsing into a single place

    Parameters
    ----------
    None

    Returns
    -------
    parser : object
        a parser object from argparse built-in python library
    '''
    parser = argparse.ArgumentParser(description=("Transcription Factor Enrichment "
                                                    "Analysis (TFEA)"))

    # The main inputs into TFEA to run the complete pipeline
    inputs = parser.add_argument_group('Main Inputs', 
                                        'Inputs required for full pipeline')
    inputs.add_argument('--output', '-o', help=("Full path to output directory. "
                        "If it exists, overwrite its contents."), dest='OUTPUT')
    inputs.add_argument('--bed1', nargs='?', help=("Bed files associated with "
                        "condition 1"), dest='BED1')
    inputs.add_argument('--bed2', nargs='?', help=("Bed files associated with "
                        "condition 2"), dest='BED2')
    inputs.add_argument('--bam1', nargs='?', help=("Sorted bam files "
                        "associated with condition 1"), 
                        dest='BAM1')
    inputs.add_argument('--bam2', nargs='?', help=("Sorted bam files "
                        "associated with condition 2"), 
                        dest='BAM2')
    inputs.add_argument('--label1', help=("An informative label for "
                        "condition 1"), dest='LABEL1')
    inputs.add_argument('--label2', help=("An informative label for "
                        "condition 2"), dest='LABEL2')
    inputs.add_argument('--config','-c', help=("A configuration file that a "
                        "user may use instead of specifying flags. Command "
                        "line flags will overwrite options within the config "
                        "file. See examples in the config_files folder."), 
                        dest='CONFIG')
    inputs.add_argument('--sbatch','-s', help=("Submits an sbatch (slurm) job. "
                        "If specified, input an e-mail address."), dest='SBATCH')
    inputs.add_argument('--test','-t', action='store_true', 
                        help="Performs unit testing", dest='TEST')

    # Processed inputs if a user desires to run TFEA from a specific point
    processed_inputs = parser.add_argument_group('Processed Inputs', 
                                        'Input options for pre-processed data')
    processed_inputs.add_argument('--combined_file', help=("A single bed file "
                                            "combining regions of interest."), 
                                            dest='COMBINED_FILE')
    processed_inputs.add_argument('--ranked_file', help=("A bed file containing each "
                                            "regions rank as the 4th column."), 
                                            dest='RANKED_FILE')
    processed_inputs.add_argument('--fasta_file', help=("A fasta file containing "
                                            "sequences to be analyzed, ranked by "
                                            "the user."), 
                                            dest='FASTA_FILE')

    # Secondary analysis inputs
    secondary_inputs = parser.add_argument_group('Secondary Analysis Inputs', 
                                                'Input options for performing '
                                                'MD-Score and Differential '
                                                'MD-Score analysis')
    secondary_inputs.add_argument('--md', help=("Switch that controls whether "
                                    "to perform MD analysis."), 
                                    action='store_true', dest='MD')
    secondary_inputs.add_argument('--mdd', help=("Switch that controls whether "
                                    "to perform differential MD analysis."), 
                                    action='store_true', dest='MDD')
    secondary_inputs.add_argument('--md_bedfile1', help=("A bed file for MD-Score "
                                    "analysis associated with condition 1."), 
                                    dest='MD_BEDFILE1')
    secondary_inputs.add_argument('--md_bedfile2', help=("A bed file for MD-Score "
                                    "analysis associated with condition 2."), 
                                    dest='MD_BEDFILE2')
    secondary_inputs.add_argument('--mdd_bedfile1', help=("A bed file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 1."), 
                                    dest='MDD_BEDFILE1')
    secondary_inputs.add_argument('--mdd_bedfile2', help=("A bed file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 2."), 
                                    dest='MDD_BEDFILE2')
    secondary_inputs.add_argument('--md_fasta1', help=("A fasta file for MD-Score "
                                    "analysis associated with condition 1."), 
                                    dest='MD_FASTA1')
    secondary_inputs.add_argument('--md_fasta2', help=("A fasta file for MD-Score "
                                    "analysis associated with condition 2."), 
                                    dest='MD_FASTA2')
    secondary_inputs.add_argument('--mdd_fasta1', help=("A fasta file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 1."), 
                                    dest='MDD_FASTA1')
    secondary_inputs.add_argument('--mdd_fasta2', help=("A fasta file for "
                                    "Differential MD-Score analysis associated "
                                    "with condition 2."), 
                                    dest='MDD_FASTA2')

    # Module switches
    module_switches = parser.add_argument_group('Module Switches', 
                                            'Switches for different modules')
    module_switches.add_argument('--combine', help=("Method for combining input "
                                    "bed files"), choices=['intersect/merge', 
                                    'merge all', 'tfit clean', 
                                    'tfit remove small', False], dest='COMBINE')
    module_switches.add_argument('--rank', help=("Method for ranking combined "
                                    "bed file"), choices=['deseq', 'fc', False], 
                                    dest='RANK')
    module_switches.add_argument('--fasta', help=("Swtich that determines "
                                    "whether converting to fasta is performed."),
                                    choices=[True, False], dest='FASTA')
    module_switches.add_argument('--scanner', help=("Method for scanning fasta "
                                    "files for motifs"), choices=['fimo', 
                                    'genome hits'], dest='SCANNER')
    module_switches.add_argument('--enrichment', help=("Method for calculating "
                                    "enrichment"), choices=['auc', 
                                    'auc_bgcorrect'], dest='ENRICHMENT')
    module_switches.add_argument('--debug', help=("Print memory usage to stderr"), 
                                    action='store_true', dest='DEBUG')
    
    # Fasta Options
    fasta_options = parser.add_argument_group('Fasta Options', ('Options for '
                                        'performing conversion of bed to fasta'))
    fasta_options.add_argument('--genomefasta', help=("Genomic fasta file"), 
                                    dest='GENOMEFASTA')

    # Scanner Options
    scanner_options = parser.add_argument_group('Scanner Options', 
                                        'Options for performing motif scanning')
    scanner_options.add_argument('--fimo_thresh', help=("P-value threshold for "
                                    "calling FIMO motif hits"), 
                                    dest='FIMO_THRESH')
    scanner_options.add_argument('--fimo_motifs', help=("Full path to a .meme "
                                    "formatted motif databse file. Some "
                                    "databases included in motif_databases "
                                    "folder."), dest='FIMO_MOTIFS')
    scanner_options.add_argument('--fimo_background', help=("Options for "
                                    "choosing mononucleotide background "
                                    "distribution to use with FIMO. "
                                    "{'largewindow', 'smallwindow', int, file}"),
                                    dest='FIMO_BACKGROUND')
    scanner_options.add_argument('--genomehits', help=("A folder containing "
                                    "bed files with pre-calculated motif hits "
                                    "to a genome. For use with 'genome hits' "
                                    "scanner option."), dest='GENOMEHITS')
    scanner_options.add_argument('--singlemotif', help=("Option to run "
                                    "analysis on a subset of motifs within "
                                    "specified motif database or genome hits. "
                                    "Can be a single motif or a comma-separated "
                                    "list of motifs."), dest='SINGLEMOTIF')

    # Enrichment Options
    enrichment_options = parser.add_argument_group('Enrichment Options', 
                                'Options for performing enrichment analysis')
    enrichment_options.add_argument('--permutations', help=("Number of "
                                        "permutations to perfrom for "
                                        "calculating p-value."), 
                                        dest='PERMUTATIONS')
    enrichment_options.add_argument('--largewindow', help=("The size (bp) of a "
                                        "large window around input regions "
                                        "that captures background."), 
                                        dest='LARGEWINDOW')
    enrichment_options.add_argument('--smallwindow', help=("The size (bp) of a "
                                        "small window arount input regions "
                                        "that captures signal."), 
                                        dest='SMALLWINDOW')
    
    # Output Options
    output_options = parser.add_argument_group('Output Options', 
                                                'Options for the output.')
    output_options.add_argument('--dpi', help=("Resolution of output figures."), 
                                    dest='DPI')
    output_options.add_argument('--padjcutoff', help=("A p-adjusted cutoff "
                                    "value that determines some plotting output."), 
                                    dest='PADJCUTOFF')
    output_options.add_argument('--pvalcutoff', help=("A p-value cutoff "
                                    "value that determines some plotting output."), 
                                    dest='PVALCUTOFF')
    output_options.add_argument('--textOnly', help="Suppress html output.", 
                                    action='store_true', dest='TEXTONLY')

    #Set default arguments and possible options or types
    # Notes: 
    # 1. 'PosixList' is a special keyword that defines a list of PosixPaths.
    #     this value is handled later in the code
    # 2. Types are casted in the order provided and handled with exceptions so
    #     it is recommended to provide the more restrictive type first.
    #
    #       CORRECT:
    #       -------
    #       'FIMO_BACKGROUND': 'largewindow', [int, str] 
    #       if FIMO_BACKGROUND = 'largewindow':
    #           int(FIMO_BACKGROUND) results in error so str(FIMO_BACKGROUND)
    #       if FIMO_BACKGROUND = 5:
    #           int(FIMO_BACKGROUND) is fine
    #
    #       INCORRECT:
    #       ---------
    #       'FIMO_BACKGROUND': 'largewindow', [str, int] 
    #       if 'FIMO_BACKGROUND' = 'largewindow':
    #           str(FIMO_BACKGROUND) is fine
    #       if FIMO_BACKGROUND = 5:
    #           str(FIMO_BACKGROUND) is fine so FIMO_BACKGROUND = '5' instead of FIMO_BACKGROUND = 5
    #
    # 3. The bool type is very permissive so if value is True/False, it is 
    #     automatically casted to bool (so specifying bool is not required)
    arg_defaults = {'OUTPUT': [False, [Path, bool]],
                    'BED1': [False, ['PosixList', bool]],
                    'BED2': [False, ['PosixList', bool]], 
                    'BAM1': [False, ['PosixList', bool]], 
                    'BAM2': [False, ['PosixList', bool]],
                    'LABEL1': [False, [str, bool]],
                    'LABEL2': [False, [str, bool]],
                    'CONFIG': [False, [Path, bool]],
                    'SBATCH': [False, [str, bool]],
                    'TEST': [False, [bool, bool]],
                    'COMBINED_FILE': [False, [Path, bool]],
                    'RANKED_FILE': [False, [Path, bool]],
                    'FASTA_FILE': [False, [Path, bool]],
                    'MD': [False, [bool]],
                    'MDD': [False, [bool]],
                    'MD_BEDFILE1': [False, [Path, bool]],
                    'MD_BEDFILE2': [False, [Path, bool]],
                    'MDD_BEDFILE1': [False, [Path, bool]],
                    'MDD_BEDFILE2': [False, [Path, bool]],
                    'MD_FASTA1': [False, [Path, bool]],
                    'MD_FASTA2': [False, [Path, bool]],
                    'MDD_FASTA1': [False, [Path, bool]],
                    'MDD_FASTA2': [False, [Path, bool]],
                    'COMBINE': ['intersect/merge', [str, bool]],
                    'RANK': ['deseq', [str, bool]],
                    'SCANNER': ['fimo', [str, bool]],
                    'ENRICHMENT': ['auc', [str]],
                    'DEBUG': [False, [bool]],
                    'GENOMEFASTA': [False, [Path, bool]],
                    'FIMO_THRESH': [0.000001, [float]], 
                    'FIMO_MOTIFS': [False, [Path, bool]],
                    'FIMO_BACKGROUND': ['largewindow', [int, str]], 
                    'SINGLEMOTIF': [False, [bool, str]], 
                    'GENOMEHITS': [False, [Path, bool]],
                    'PERMUTATIONS': [1000, [int]], 
                    'LARGEWINDOW': [1500, [int]], 
                    'SMALLWINDOW': [150, [int]], 
                    'DPI': [None, [str, int]], 
                    'PADJCUTOFF': [0.001, [float]], 
                    'PVALCUTOFF': [0.01, [float]], 
                    'TEXTONLY': [False, [bool]]}

    #Save default arguments in config
    import config
    config.arg_defaults = arg_defaults

    return parser

#==============================================================================
def make_out_directories(create=False):
    '''Creates output directories in a user-specified location where all TFEA 
        outputs will go.

    Parameters
    ----------
    create : boolean
        determines whether output folders will be created or not (default: 
        False)
        
    config_dict : dict
        a configparser object that contains variables within the config file

    Returns
    -------
    output : string 
        full path to the parent output directory

    tempdir : string
        full path to the temporary directory where files are stored

    figuredir : string
        full path to the directory containing figures and plots
        
    e_and_o : string
        full path to the directory that stores stdout and stderr files
    '''
    import config
    #Output directory
    output = config.vars.OUTPUT

    #Make parent output directory
    if create:
        output.mkdir(exist_ok=True)

    #Temporary files will go in this directory
    tempdir = output / 'temp_files'
    if create:
        tempdir.mkdir(exist_ok=True)

    #Error and out files will go in this directory
    e_and_o = output / 'e_and_o'
    if create:
        e_and_o.mkdir(exist_ok=True)

    #Directory where plots used in html file will be stored.
    figuredir = output / 'plots'
    if create:
        figuredir.mkdir(exist_ok=True)

    config.vars.TEMPDIR = tempdir
    config.vars.FIGUREDIR = figuredir
    config.vars.E_AND_O = e_and_o

#==============================================================================
def parse_config(configfile=None):
    '''This function parses a config file converting it into a configparser
        object and converting that to a python dictionary

    Parameters
    ----------
    configfile : str
        full path to a config file (.ini)
    
    Returns
    -------
    config_dict : dict
        a dict where keys are variables and values are user inputs. Keys are
        converted to uppercase
    '''
    import configparser
    config_object = configparser.ConfigParser(
                            interpolation=configparser.ExtendedInterpolation())
    config_object.read(configfile)
    config_dict = dict()
    for key in config_object:
        for item in config_object[key]:
            config_dict[item.upper()] = eval(config_object[key][item])

    return config_dict

#==============================================================================
def verify_arguments(parser=None):
    '''Verifies that all necessary variables are present within the inputted 
        config file and that they are the correct variable types.
    
    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    InputError
        When missing or conflicting input
    '''
    import config
    #If a config file is specified, parse its arguments
    configfile_dict = dict()
    configfile = parser.parse_args().CONFIG
    if configfile != None:
        configfile_dict = parse_config(configfile=configfile)

    #Overwrite config file args with command line flag args
    config_dict = dict()
    parser_dict = vars(parser.parse_args())
    arg_defaults = config.arg_defaults
    for key in arg_defaults:
        #Decide which value to use
        if key in configfile_dict and parser_dict[key] == None:
            value = configfile_dict[key] #User did not specify a flag but specified variable in config file
        elif parser_dict[key] != None:
            value = parser_dict[key] #User specified a flag (overwrites config file)
        else:
            value = arg_defaults[key][0] #Nothing specified, use default argument

        #Cast to correct type
        if value in (True, 'True', False, 'False'): #Detect bools explicitly
            config_dict[key] = bool(value)
        elif 'PosixList' in arg_defaults[key][1]: #Handle special 'PosixList' type
            types = arg_defaults[key][1]
            for default_type in types:
                try:
                    value = [Path(x) for x in value]
                    break
                except TypeError:
                    try:
                        value = default_type(value)
                        break
                    except TypeError:
                        pass
            config_dict[key] = value
        else: # For all other types
            types = arg_defaults[key][1]
            for default_type in types:
                try:
                    value = default_type(value)
                    break
                except (TypeError, ValueError):
                    pass
            config_dict[key] = value

    #Write variables to config
    config.vars = Bunch(config_dict)
    config.vars.COMBINEtime = 0
    config.vars.COUNTtime = 0
    config.vars.RANKtime = 0
    config.vars.FASTAtime = 0
    config.vars.SCANNERtime = 0
    config.vars.ENRICHMENTtime = 0
    config.vars.OUTPUTtime = 0

    #Verify combine module
    if not config.vars.COMBINE:
        if not config.vars.COMBINED_FILE and config.vars.RANK:
            raise exceptions.InputError('COMBINE module switched off but RANK module switched on without a COMBINED_FILE')
        if config.vars.MD:
            if not config.vars.MD_BEDFILE1:
                raise exceptions.InputError('COMBINE module switched off but MD module switched on without MD_BEDFILE1')
            if not config.vars.MD_BEDFILE2:
                raise exceptions.InputError('COMBINE module switched off but MD module switched on without MD_BEDFILE2')
    else:
        if not config.vars.BED1:
            raise exceptions.InputError('COMBINE module switched on but BED1 not specified')
        if not config.vars.BED2:
            raise exceptions.InputError('COMBINE module switched on but BED2 not specified')

    #Verify rank module
    if not config.vars.RANK:
        if not config.vars.RANKED_FILE and not config.vars.FASTA_FILE and config.vars.scanner == 'fimo':
            raise exceptions.InputError('SCANNER module set to "fimo" but RANK module switched off without RANKED_FILE or FASTA_FILE')
        if config.vars.MDD:
            if not config.vars.MDD_BEDFILE1:
                raise exceptions.InputError('RANK module switched off but MDD module switched on without MDD_BEDFILE1')
            if not config.vars.MDD_BEDFILE2:
                raise exceptions.InputError('RANK module switched off but MDD module switched on without MDD_BEDFILE2')
    else:
        if not config.vars.BED1:
            raise exceptions.InputError('RANK module switched on but BED1 not specified')
        if not config.vars.BED2:
            raise exceptions.InputError('RANK module switched on but BED2 not specified')
        if not config.vars.BAM1:
            raise exceptions.InputError('RANK module switched on but BAM1 not specified')
        if not config.vars.BAM2:
            raise exceptions.InputError('RANK module switched on but BAM2 not specified')
        if not config.vars.LABEL1:
            raise exceptions.InputError('RANK module switched on but LABEL1 not specified')
        if not config.vars.LABEL2:
            raise exceptions.InputError('RANK module switched on but LABEL2 not specified')
        

    #Verify scanner module
    if not config.vars.GENOMEHITS and config.vars.SCANNER == 'genome hits':
        raise exceptions.InputError('SCANNER set to "genome hits" without specifying GENOMEHITS')

    print("User arguments verified, all required inputs present and not conflicting.", 
            file=sys.stderr)

#==============================================================================
def create_directories(srcdirectory=None):
    import config
    if config.vars.SBATCH == False: #No sbatch flag
        make_out_directories(create=True)
        write_rerun(args=sys.argv, outputdir=config.vars.OUTPUT)
    elif str(config.vars.SBATCH ) == 'SUBMITTED': #Internal flag
        make_out_directories(create=False)
    else: #--sbatch specified
        make_out_directories(create=True)
        write_rerun(args=sys.argv, outputdir=config.vars.OUTPUT)
        scriptdir = srcdirectory.parent / 'cluster_scripts'
        script = scriptdir / 'run_main.sbatch'
        email = str(config.vars.SBATCH)
        error_file = config.vars.E_AND_O / 'TFEA.err'
        args = sys.argv
        args[args.index('--sbatch')+1] = 'SUBMITTED'
        try:
            sbatch_out = subprocess.run(["sbatch", 
                                "--error=" + (config.vars.E_AND_O / "%x.err").as_posix(), 
                                "--output=" + (config.vars.E_AND_O / "%x.out").as_posix(), 
                                "--mail-user=" + email, 
                                "--export=cmd=" + ' '.join(args), 
                                script], stderr=subprocess.PIPE, 
                                stdout=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())
        
        (config.vars.TEMPDIR / 'jobid.txt').write_text(sbatch_out.stdout.decode().split()[-1])
        print(("TFEA has been submitted using an sbatch script. \nIt can be "
                "monitored using:\ntail -f " + error_file.as_posix()))
        sys.exit()

#==============================================================================
def write_rerun(args=None, outputdir=None):
    with open(outputdir / 'rerun.sh', 'w') as outfile:
        outfile.write(' '.join(args))
    try:
        configfile = args[args.index('--config')+1]
        try:
            shutil.copy2(configfile, outputdir) #Copy config file to output
        except shutil.SameFileError:
            pass
    except ValueError:
        pass

