#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module contains several basic test cases to determine whether TFEA
    is installed properly and working.
'''

#==============================================================================
__author__ = 'Jonathan D. Rubin and Rutendo F. Sigauke'
__credits__ = ['Jonathan D. Rubin', 'Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Jonathan D. Rubin'
__email__ = 'Jonathan.Rubin@colorado.edu'

#Imports
#==============================================================================
import os
import sys
import unittest
import subprocess
import shutil
from pathlib import Path

#Tests
#==============================================================================
def execute(cmd):
    '''Function to execute a command using subprocess and have the output
        printed in real time. Taken from a StackOverflow question:
        (https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running#)
    '''
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

class TestMain(unittest.TestCase):
    def setUp(self):
        self.srcdir = Path(__file__).parent
        self.testdir = self.srcdir / 'test_files'
        self.motifs_meme = self.testdir / 'test_database.meme'
        self.include_motifs = 'SP2_HUMAN.H11MO.0.A'
        self.random_fasta = self.testdir / 'Allen2014_simulated_10.fa'
        self.sequence_n = '5'
        self.smallwindow = '150'
        self.largewindow = '1500'
        self.distance_mu = '0'
        self.distance_sigma = 'uniform'
        self.rank_range = '0-5'
        self.motif_number = '5'
        self.seed = '0'

    def test_simulate(self):
        simulate_src_path = Path(__file__).absolute().parent.parent / 'simulate'
        command = ['nice', '-n', '19',
                    'python3', simulate_src_path, 
                    '--output', self.testdir / 'test_output' / 'simulate_test.fa', 
                    '--motifs', self.motifs_meme, 
                    '--include_motifs', self.include_motifs, 
                    '--random_fasta', self.random_fasta, 
                    '--sequence_n', self.sequence_n, 
                    '--smallwindow', self.smallwindow, 
                    '--largewindow', self.largewindow,
                    '--distance_mu', self.distance_mu,
                    '--distance_sigma', self.distance_sigma,  
                    '--rank_range', self.rank_range,
                    '--motif_number', self.motif_number,
                    '--seed', self.seed]

        print('\n============================================', 
                file=sys.stderr, flush=True)
        for output in execute(command):
            try:
                print('\t', output.decode(), end="", flush=True)
            except:
                print('\t', output, end="", flush=True)

        self.assertTrue((self.testdir / 'test_output' / 'simulate_test.fa').exists())
        # self.assertTrue((self.testdir / 'test_output' / 'results.html').exists())
        # self.assertTrue((self.testdir / 'test_output' / 'summary.html').exists())
        # self.assertTrue((self.testdir / 'test_rep1' / 'md_results.txt').exists())
        # self.assertTrue((self.testdir / 'test_rep1' / 'mdd_results.txt').exists())

if __name__ == '__main__':
    import sys
    from pathlib import Path
    # Add TFEA srcdirectory into path
    srcdirectory = Path(__file__).absolute().parent.parent
    sys.path.insert(0, srcdirectory)
    unittest.main(verbosity=2, failfast=True)