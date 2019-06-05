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
import warnings
import tracemalloc
import subprocess
import shutil
from pathlib import Path

import numpy as np

from .context import exceptions

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
        self.srcdir = Path(__file__).parent.parent.parent
        self.testdir = self.srcdir / 'test_files'
        self.bed1 = [self.testdir / 'SRR1105736.tfit_bidirs.chr22.bed', 
                    self.testdir / 'SRR1105737.tfit_bidirs.chr22.bed']
        self.bed2 = [self.testdir / 'SRR1105738.tfit_bidirs.chr22.bed', 
                    self.testdir / 'SRR1105739.tfit_bidirs.chr22.bed']
        self.smallwindow = 150
        self.largewindow = 1500
        self.fimo_thresh = 0.000001
        self.bam1 = [self.testdir / 'SRR1105736.sorted.chr22.subsample.bam', 
                    self.testdir / 'SRR1105737.sorted.chr22.subsample.bam']
        self.bam2 = [self.testdir / 'SRR1105738.sorted.chr22.subsample.bam', 
                    self.testdir / 'SRR1105739.sorted.chr22.subsample.bam']
        self.label1 = 'condition1'
        self.label2 = 'condition2'
        self.genomefasta = self.testdir / 'chr22.fa'
        self.ranked_file = self.testdir / 'ranked_file.bed'

        self.count_file = self.testdir / 'count_file.header.bed'
        self.fimo_motifs = self.testdir / 'test_database.meme'

        touch_command = "touch " + str(self.testdir / '*.bai')
        subprocess.run(touch_command, check=True, shell=True)
    
    def test_full(self):
        shutil.rmtree(self.testdir / 'test_output', ignore_errors=True)
        TFEA_path = Path(__file__).absolute().parent.parent
        command = ['nice', '-n', '19',
                    'python3', TFEA_path, 
                    '--output', self.testdir / 'test_output', 
                    '--bed1', self.bed1[0], 
                    '--bed2', self.bed2[0], 
                    '--bam1', self.bam1[0], 
                    '--bam2', self.bam2[0], 
                    '--label1', self.label1, 
                    '--label2', self.label2,
                    '--genomefasta', self.genomefasta,
                    '--fimo_motifs', self.fimo_motifs,  
                    '--motif_annotation', self.testdir / 'test_motif_annotation.bed',
                    '--output_type', 'html', 
                    '--plotall',
                    '--combine', 'merge all',
                    '--md', '--mdd',
                    '--debug']
        # try:
        print('\n============================================')
        for output in execute(command):
            print('\t', output.decode(), end="", flush=True)
            # subprocess.check_output(command, stderr=subprocess.PIPE)
        # except subprocess.CalledProcessError as e:
        #     raise exceptions.SubprocessError(e.stderr.decode())

        self.assertTrue((self.testdir / 'test_output' / 'results.txt').exists())
        self.assertTrue((self.testdir / 'test_output' / 'results.html').exists())
        # self.assertTrue((self.testdir / 'test_output' / 'summary.html').exists())
        # self.assertTrue((self.testdir / 'test_rep1' / 'md_results.txt').exists())
        # self.assertTrue((self.testdir / 'test_rep1' / 'mdd_results.txt').exists())

    # def test_2rep(self):
    #     shutil.rmtree(self.testdir / 'test_rep2', ignore_errors=True)
    #     TFEA_path = Path(__file__).absolute().parent.parent
    #     command = ['nice', '-n', '19',
    #                 'python3', TFEA_path, '--output', self.testdir / 'test_rep2', 
    #                 '--bed1'] + self.bed1 + ['--bed2'] +  self.bed2 + ['--bam1'] + self.bam1 + ['--bam2'] + self.bam2 + ['--label1', self.label1, 
    #                 '--label2', self.label2,
    #                 '--genomefasta', self.genomefasta,
    #                 '--fimo_motifs', self.fimo_motifs,  
    #                 '--output_type', 'html', 
    #                 '--combine', 'merge all',
    #                 '--plotall',
    #                 '--debug',
    #                 '--md', '--mdd']
    #     # try:
    #     print()
    #     for output in execute(command):
    #         print('\t', output, end="")
    #         # subprocess.check_output(command, stderr=subprocess.PIPE)
    #     # except subprocess.CalledProcessError as e:
    #     #     raise exceptions.SubprocessError(e.stderr.decode())

    #     self.assertTrue((self.testdir / 'test_rep2' / 'results.txt').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'results.html').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'summary.html').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'md_results.txt').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'mdd_results.txt').exists())

if __name__ == '__main__':
    unittest.main()