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
import unittest
import warnings
import tracemalloc
from pathlib import Path

from .context import config
from .context import combine
from .context import count
from .context import rank
from .context import fasta
from .context import scanner
from .context import enrichment
from .context import multiprocess

#Tests
#==============================================================================
class TestMain(unittest.TestCase):
    def setUp(self):
        self.srcdir = Path(__file__).parent.parent.parent
        self.testdir = self.srcdir / 'test_files'
        self.bed1 = [self.testdir / 'bed1_rep1.bed', 
                    self.testdir / 'bed1_rep2.bed']
        self.bed2 = [self.testdir / 'bed2_rep1.bed', 
                    self.testdir / 'bed2_rep2.bed']
        self.smallwindow = 150
        self.largewindow = 1500
        self.bedfile = self.testdir / 'bed1_rep1.bed'
        self.bam1 = [self.testdir / 'bed1_rep1.bam', 
                    self.testdir / 'bed1_rep2.bam']
        self.bam2 = [self.testdir / 'bed2_rep1.bam', 
                    self.testdir / 'bed2_rep2.bam']
        self.label1 = 'condition1'
        self.label2 = 'condition2'
        self.count_file = self.testdir / 'test_count_file.bed'
        self.ranked_file = self.testdir / 'test_bed.sorted.bed'
        self.genomefasta = self.testdir / 'test_genome.fa'
        self.fasta_file = self.testdir / 'test_fasta_file.fa'
        self.fimo_motifs = self.testdir / 'test_database.meme'
        self.combined_file = self.testdir / 'test_combined_file.bed'
    
    def test_mergeall(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bedfile = combine.main(bed1=self.bed1, bed2=self.bed2, 
                                    method='merge all', tempdir=self.testdir, 
                                    md=False, largewindow=self.largewindow, 
                                    scanner='fimo')
            self.assertEqual(len(config.vars.COMBINED_FILE.read_text().strip('\n').split('\n')), 11)

    def test_intersectall(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            combine.main(bed1=[self.bed1[0]], bed2=[self.bed2[0]], 
                                    method='intersect all', tempdir=self.testdir, 
                                    md=False, largewindow=self.largewindow, 
                                    scanner='fimo')
            self.assertEqual(len(config.vars.COMBINED_FILE.read_text().strip('\n').split('\n')), 1)

    def test_intersectmerge(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bedfile = combine.main(bed1=self.bed1, bed2=self.bed2, 
                                    method='intersect/merge', tempdir=self.testdir, 
                                    md=False, largewindow=self.largewindow, 
                                    scanner='fimo')
            self.assertEqual(len(config.vars.COMBINED_FILE.read_text().strip('\n').split('\n')), 9)

    # def test_count(self):
    #     count_file = count.main(bedfile=self.bedfile, bam1=self.bam1, 
    #                         bam2=self.bam2, 
    #                         tempdir=self.testdir, label1=self.label1, 
    #                         label2=self.label2)
    #     self.assertEqual(len(count_file.read_text().strip('\n').split('\n')), 7)

    # def test_rank(self):
    #     ranked_file = rank.main(combined_file=self.combined_file, rank='deseq', scanner='fimo', 
    #         bam1=self.bam1, bam2=self.bam2, tempdir=self.testdir, 
    #         label1=self.label1, label2=self.label2, 
    #         largewindow=self.largewindow)
    #     self.assertEqual(len(ranked_file.read_text().strip('\n').split('\n')), 11)

    # def test_fasta(self):
    #     fasta_file = fasta.main(ranked_file=self.ranked_file, 
    #         genomefasta=self.genomefasta, tempdir=self.testdir)
    #     self.assertEqual(len(fasta_file.read_text().strip('\n').split('\n')), 8)

    def test_scanner(self):
        motif_distances = scanner.main(fasta_file=self.fasta_file, 
                                        ranked_file=self.ranked_file, 
                                        scanner='fimo', md=False, 
                                        largewindow=self.largewindow,
                                        smallwindow=self.smallwindow, 
                                        genomehits=None, 
                                        fimo_background=None, 
                                        genomefasta=self.genomefasta, 
                                        tempdir=self.testdir, 
                                        fimo_motifs=self.fimo_motifs, 
                                        singlemotif=False, 
                                        fimo_thresh='1e-6', debug=False)
        self.assertEqual(len(config.vars.MOTIF_DISTANCES[0]), 5)

    def tearDown(self):
        [path.unlink() for path in self.testdir.iterdir() if 'bed1' not in path.name and 'bed2' not in path.name and 'test_' not in path.name]

if __name__ == '__main__':
    unittest.main()