#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module contains scripts to perform parallel processing using python's
    built-in multiprocessing module
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import multiprocessing as mp
import subprocess
import resource
import psutil
import sys
from pathlib import Path
from contextlib import closing

from TFEA import config

#Main Script
#==============================================================================
def main(function=None, args=None, kwargs=None, debug=False, jobid=None, cpus=None):
    '''This is the main script of the multiprocessing module. It is written
        to simplify code in other modules. It performs parallel processing
        using python's built-in multiprocessing module on a function given
        keyword arguments shared between threads and arguments given to each
        thread

    Parameters
    ----------
    function : function object
        A python function object to process some arguments
    args : list
        A list of arguments to be input into given function
    kwargs : dict
        A dictionary of keyword arguments to input into given function
    debug : boolean
        Whether to print memory usage information of running processes

    Returns
    -------
    results : list
        A list of the result of processing given args and kwargs with given 
        function
    '''
    # #Get max number of cpus, max out at 64
    # cpus = mp.cpu_count()
    # if cpus > 64:
    #     cpus = 64
    # if debug:
    #     print('CPUs:', cpus)

    #Method 1
    # #with closing(mp.Pool(processes=cpus)) as p:
    # with mp.get_context("spawn").Pool(maxtasksperchild=10, processes=cpus) as p:
    #     results = p.starmap(helper, [(function, arg, kwargs, debug, i, jobid) for i, arg in enumerate(args, 1)])
    #     p.close()
    #     p.join()

    #Method 2 
    with mp.Pool(cpus) as p:
        results = list()
        print(f'\t Completed: 0/{len(args)} ', end=' ', file=sys.stderr)
        # print_in_place(f'\t Completed: 0/{len(args)} ', file=sys.stderr)
        for i, x in enumerate(p.imap_unordered(helper_single, [(function, arg, kwargs, debug) for arg in args]), 1):
            # sys.stderr.write("\033[K")
            # sys.stderr.write(f'\r\t Completed: {i}/{len(args)} ')
            # sys.stderr.flush()
            print(f'\r\t Completed: {i}/{len(args)} ', end=' ', flush=True, file=sys.stderr)
            # print_in_place(f'\t Completed: {i}/{len(args)} ', file=sys.stderr)
            if debug:
                current_mem_usage(jobid, end=' ')
                # current_mem_usage(jobid, in_place=True)
            results.append(x)
        p.close()
        p.join()
    print('', file=sys.stderr)

    #Method 3
    # with mp.get_context("spawn").Pool(maxtasksperchild=10, processes=cpus) as p:
    #     results = p.imap_unordered(helper_single, [(function, arg, kwargs, debug, i, jobid) for i, arg in enumerate(args, 1)])
    #     print('\n', end=' ', file=sys.stderr)
    #     for i in range(len(args)):
    #         results.next()
    #         sys.stderr.write("\033[K")
    #         sys.stderr.write(f'\r\tCompleted: {i}/{len(args)} ')
    #         sys.stderr.flush()
    #         if debug:
    #             current_mem_usage(jobid, end=' ', flush=True)
    #     p.close()
    #     p.join()
    #     p.terminate()

    
        
    # for i, arg in enumerate(args):
    #     results.append(p.apply_async(function, (arg,), kwargs).get())
    #     if debug:
    #         print("Process number: " + str(i), file=sys.stderr)
    #         current_mem_usage()
    # p.close()
    # p.join()
    # p.terminate()

    return results

#Functions
#==============================================================================
def helper(function, arg, kwargs, debug, i, jobid):
    '''This function serves to unpack keyword arguments provided to this module
        and feed it in appropriately to the desired function

    Parameters
    ----------
    function : function object
        A python function object to process some arguments
    arg : object
        A positional argument to pass to function
    kwargs : dict
        keyword arguments to pass to function
    debug : boolean
        whether to print debug statements
    i : int
        the process number
    jobid : int
        the slurm pid

    Returns
    -------
    result : object
        The result of running the function with given argument and keyword
        arguments.
    '''
    if debug:
        print(f"Process #{i}, pid:{os.getpid()} running {function.__name__}",
                file=sys.stderr)
        current_mem_usage(jobid)
    result = function(arg, **kwargs)
    if debug:
        # print("Process", os.getpid(), 
        #         "finished", function.__name__, file=sys.stderr)
        print(f"Process #{i}, pid:{os.getpid()} finished {function.__name__}",
                file=sys.stderr)

    return result

#==============================================================================
def helper_single(args):
    '''This function serves to unpack keyword arguments provided to this module
        and feed it in appropriately to the desired function. Different than
        helper in that it takes a single argument.

    Parameters
    ----------
    args : tuple
        Contains all necessary arguments to run helper_single

    Returns
    -------
    result : object
        The result of running the function with given argument and keyword
        arguments.
    '''
    function, arg, kwargs, debug = args
    # if debug:
    #     print(f"Process #{i}, pid:{os.getpid()} running {function.__name__}",
    #             file=sys.stderr)
    #     current_mem_usage(jobid)
    result = function(arg, **kwargs)
    # if debug:
    #     # print("Process", os.getpid(), 
    #     #         "finished", function.__name__, file=sys.stderr)
    #     print(f"Process #{i}, pid:{os.getpid()} finished {function.__name__}",
    #             file=sys.stderr)

    return result

#==============================================================================
def current_mem_usage(jobid, **kwargs):
    '''Prints current memory usage. Adapted from scripts by Chris Slocum and 
        Fred Cirera. Additional kwargs are passed to print function
    '''
    suffix = 'B'

    #Using sstat
    # num = subprocess.check_output(["sstat", "--format=maxRSS", "-j", jobid+'.batch']).decode('UTF-8').split('\n')[-2].strip()
    # unit = num[-1]

    #Using psutil
    num = psutil.Process(pid=os.getpid()).memory_info().rss
    # num = psutil.virtual_memory().active
    unit = 'B'
    print(f'(CPU: {psutil.cpu_percent()}%)', file=sys.stderr, **kwargs)

    #Using resource
    # num = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # unit = 'B'

    if unit == 'B':
        byte = float(str(num).strip('B'))
    if unit == 'K':
        byte = float(str(num).strip('K'))*1024.0
    elif unit == 'M':
        byte = float(str(num).strip('M'))*1024.0*1024.0
    elif unit == 'G':
        byte = float(str(num).strip('G'))*1024.0*1024.0*1024.0
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(byte) < 1024.0:
            print("(Memory Usage: %3.1f%s%s)" % (byte, unit, suffix), file=sys.stderr, **kwargs)
            return
        byte /= 1024.0
    else:
        print("(Memory Usage: %.1f%s%s)" % (byte, 'Y', suffix), file=sys.stderr, **kwargs)
    return

#==============================================================================
def limit_cpu():
    "is called at every process start"
    p = psutil.Process(os.getpid())
    p.nice(value=19)