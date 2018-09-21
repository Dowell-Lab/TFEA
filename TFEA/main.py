def run():
    #==============================================================================
    #MAIN IMPORTS
    #==============================================================================
    import os
    import sys
    import time
    import argparse
    import configparser
    #==============================================================================
    #ARGUMENT PARSING
    #==============================================================================
    #TFEA source directory
    srcdirectory = os.path.dirname(os.path.realpath(__file__))

    #argparse to add arguments to this python package
    parser = argparse.ArgumentParser(description='Transcription Factor Enrichment \
                                        Analysis (TFEA) takes as input a \
                                        configuration file (.ini) and outputs \
                                        a folder containing TFEA results.',
                                    usage='TFEA --config CONFIG.ini [--sbatch \
                                        email@address.com]')

    parser.add_argument('--config','-c',metavar='',help='REQUIRED. A \
                            configuration file containing .ini suffix \
                            (ex. config.ini). See example in the examples folder.')

    parser.add_argument('--sbatch','-s',default=False,metavar='',help='OPTIONAL. \
                            Submits an sbatch job. If specified, input an e-mail \
                            address.')

    #Display help message when no args are passed.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    #==============================================================================
    #SECONDARY IMPORTS
    #==============================================================================
    from multiprocessing import Pool
    import multiprocessing as mp
    import independent_functions
    import dependent_functions
    #==============================================================================
    #MAIN SCRIPT
    #==============================================================================
    #If user provided arguments, then parse them
    sbatch = parser.parse_args().sbatch
    configfile = parser.parse_args().config
    config_object = configparser.ConfigParser(
                                interpolation=configparser.ExtendedInterpolation())
    config_object.read(configfile)

    #If user specifies the --sbatch flag, then we first create the output 
    #directories then run the sbatch script with the 'SUBMITTED' command submitted 
    #to the --sbatch flag so we know not to remake output directories. If --sbatch 
    #flag not specified, simply make output directories and continue.
    if sbatch == False:
        args = independent_functions.make_out_directories(dirs=True, 
                                                            config=config_object)
        output, tempdir, figuredir, e_and_o = args
    elif str(sbatch) == 'SUBMITTED':
        args = independent_functions.make_out_directories(dirs=False, 
                                                            config=config_object)
        output, tempdir, figuredir, e_and_o = args
    else:
        args = independent_functions.make_out_directories(dirs=True, 
                                                            config=config_object)
        output, tempdir, figuredir, e_and_o = args
        scriptdir = os.path.join(os.path.dirname(srcdirectory), 'scripts')
        script = os.path.join(scriptdir, 'run_main.sbatch')
        email = str(sbatch)
        os.system("sbatch --error=" + e_and_o + "/%x.err --output=" + e_and_o 
                    + "/%x.out --mail-user=" + email + " --export=cmd='" 
                    + srcdirectory+" --config " + configfile 
                    + " --sbatch SUBMITTED' " + script)

        sys.exit(("TFEA has been submitted using an sbatch script, use qstat to "
                "check its progress."))


    #Run the config_parser script which will create variables for all folders and 
    #paths to use throughout TFEA
    print "srcdirectory: ", srcdirectory
    independent_functions.parse_config(srcdirectory=srcdirectory, 
                                        config_object=config_object,
                                        output=output,tempdir=tempdir,
                                        figuredir=figuredir)

    #Verify config file to make sure user has inputted all necessary variables
    # config_dict = independent_functions.verify_config_object(config=config)
    independent_functions.verify_config_file()

    #Import config file once it's created
    import config
    import sys

    print "config temp: ", config.TEMPDIR
    config.TEMPDIR = tempdir
    print "config temp new: ", config.TEMPDIR
    print "main temp: ", tempdir

    sys.exit(1)