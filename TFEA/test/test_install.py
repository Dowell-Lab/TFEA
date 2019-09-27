#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This module is used to determine whether all required programs are installed
    and command-line runnable.
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
import unittest
import subprocess
import warnings
from pathlib import Path

#Tests
#==============================================================================
class TestMain(unittest.TestCase):
    def setUp(self):
        self.srcdir = Path(__file__).parent
        self.testdir = self.srcdir / 'test_files'
        # self.testdir = Path(__file__).parent

    def test_python_packages(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            import matplotlib
            matplotlib.use('Agg')
            import numpy
            import scipy
            import pybedtools
            import psutil
            import HTSeq

    def test_bedtools(self):
        command = ["bedtools", "--help"]
        output = subprocess.run(command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('command not found' in output.stdout.decode())
        except AssertionError:
            raise AssertionError('Bedtools not installed')

    def test_meme(self):
        command = ["meme", "-version"]
        output = subprocess.run(command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('command not found' in output.stdout.decode())
        except AssertionError:
            raise AssertionError('MEME not installed')
    
    def test_fimo(self):
        command = ["fimo", "--version"]

        output = subprocess.run(command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('command not found' in output.stdout.decode())
        except AssertionError:
            raise AssertionError('FIMO not runnable, check MEME installation')
    
    def test_fasta_markov(self):
        command = ["fasta-get-markov", "--help"]
        
        output = subprocess.run(command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('command not found' in output.stdout.decode())
        except AssertionError:
            raise AssertionError('fasta-get-markov not runnable, check MEME installation')

    def test_meme2images(self):
        command = ["meme2images", "--help"]

        output = subprocess.run(command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('command not found' in output.stdout.decode())
        except AssertionError:
            raise AssertionError(('meme2images not runnable. Motif logos will '
                                    'not be displayed in html output. Check '
                                    'MEME installation to fix.'))

    def test_imagemagick(self):
        command = ["identify", "-version"]

        output = subprocess.run(command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('command not found' in output.stdout.decode())
        except AssertionError:
            raise AssertionError(('ImageMagick not installed. Motif logos will '
                                    'not be displayed in html output'))

    def test_deseq(self):
        Rscript = Path(self.testdir / 'deseq_test.R')
        Rscript.write_text('''library(DESeq)''')
        R_command = ["Rscript", "--vanilla", Rscript]

        output = subprocess.run(R_command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('there is no package called' in output.stdout.decode())
        except AssertionError:
            raise AssertionError("DE-Seq not installed")

    def test_deseq2(self):
        Rscript = Path(self.testdir / 'deseq_test.R')
        Rscript.write_text('''library(DESeq2)''')
        R_command = ["Rscript", "--vanilla", Rscript]
        output = subprocess.run(R_command, stdout=subprocess.PIPE, 
                                stderr=subprocess.STDOUT)
        try:
            self.assertFalse('there is no package called' in output.stdout.decode())
        except AssertionError:
            raise AssertionError("DE-Seq2 not installed")

if __name__ == '__main__':
    unittest.main(verbosity=2)