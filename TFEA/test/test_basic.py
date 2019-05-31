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
from pathlib import Path

import numpy as np

from .context import config
from .context import combine
from .context import rank
from .context import scanner
from .context import enrichment
from .context import multiprocess
from .context import exceptions

#Tests
#==============================================================================
class TestModules(unittest.TestCase):
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
        self.fasta_file = self.testdir / 'test_fasta_file.fa'
        self.fimo_motifs = self.testdir / 'test_database.meme'
        self.combined_file = self.testdir / 'combined_file.mergeall.bed'
        self.motif_distances = [['motif'] + [0 for i in range(50)] + [500 for i in range(50)] + ['.' for i in range(50)]]
        self.pvals = list(np.linspace(0, 1, num=75)) + list(np.linspace(1, 0, num=75))
        self.fcs = list(np.linspace(2, 1, num=75)) + list(np.linspace(1, 0.1, num=75))
        self.profile = list(np.linspace(0, 500, 3000))
        self.metaprofile = dict({'q1posprofile1': self.profile, 
                        'q1negprofile1': self.profile,
                        'q1posprofile2': self.profile, 
                        'q1negprofile2': self.profile, 
                        'q2posprofile1': self.profile, 
                        'q2negprofile1': self.profile, 
                        'q2posprofile2': self.profile, 
                        'q2negprofile2': self.profile,
                        'q3posprofile1': self.profile, 
                        'q3negprofile1': self.profile, 
                        'q3posprofile2': self.profile, 
                        'q3negprofile2': self.profile,
                        'q4posprofile1': self.profile, 
                        'q4negprofile1': self.profile, 
                        'q4posprofile2': self.profile, 
                        'q4negprofile2': self.profile})

        # config.vars.BED1 = self.bed1
        # config.vars.BED2 = self.bed2
        # config.vars.TEMPDIR = self.testdir
        # config.vars.LARGEWINDOW = self.largewindow
    
    def test_mergeall(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            combine.main(use_config=False, bed1=self.bed1, bed2=self.bed2, 
                                    method='merge all', tempdir=self.testdir, 
                                    md=False, largewindow=self.largewindow, 
                                    scanner='fimo')
            # combine.main(use_config=True)
            self.assertEqual(len((self.testdir / "combined_file.mergeall.bed").read_text().strip('\n').split('\n')), 371)

    def test_intersectall(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            combine.main(use_config=False, bed1=self.bed1, bed2=self.bed2, 
                                    method='intersect all', tempdir=self.testdir, 
                                    md=False, largewindow=self.largewindow, 
                                    scanner='fimo')
            self.assertEqual(len((self.testdir / 'combined_file.intersectall.bed').read_text().strip('\n').split('\n')), 51)

    def test_intersectmerge(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            combine.main(use_config=False, bed1=self.bed1, bed2=self.bed2, 
                                    method='intersect/merge', tempdir=self.testdir, 
                                    md=False, largewindow=self.largewindow, 
                                    scanner='fimo')
            self.assertEqual(len((self.testdir / 'combined_file.intermerge.bed').read_text().strip('\n').split('\n')), 110)

    def test_rank(self):
        rank.main(use_config=False, combined_file=self.combined_file, 
                    rank='deseq', scanner='fimo', bam1=self.bam1, 
                    bam2=self.bam2, tempdir=self.testdir, label1=self.label1, 
                    label2=self.label2, largewindow=self.largewindow)
        self.assertEqual(len(self.ranked_file.read_text().strip('\n').split('\n')), 372)

    def test_scanner(self):
        motif_distances, _, _, _, _ = scanner.main(use_config=False, 
                                                ranked_file=self.ranked_file, 
                                                scanner='fimo', md=False, 
                                                largewindow=self.largewindow, 
                                                smallwindow=self.smallwindow, 
                                                fimo_background='largewindow', 
                                                genomefasta=self.genomefasta, 
                                                tempdir=self.testdir, 
                                                fimo_motifs=self.fimo_motifs, 
                                                singlemotif=False, 
                                                fimo_thresh=self.fimo_thresh,
                                                debug=False, mdd=False, cpus=1)
        self.assertEqual(len(motif_distances[0]), len(self.ranked_file.read_text().strip('\n').split('\n')))
    
    def test_enrichment(self):
        results, _, _ = enrichment.main(use_config=False, motif_distances=self.motif_distances,
            enrichment='auc', output_type='html', permutations=1000, debug=False, 
            largewindow=self.largewindow, smallwindow=self.smallwindow, md=False, 
            mdd=False, figuredir=self.testdir, fimo_motifs=self.fimo_motifs,
            cpus=1, p_cutoff=0.01, pvals=self.pvals, fcs=self.fcs, 
            label1=self.label1, label2=self.label2, dpi=100, meta_profile_dict=self.metaprofile)
        self.assertEqual(results[0][2], 100)

    # def test_multiprocess(self):
    #     def test_function(size):
    #         results = list()
    #         for i in range(size):
    #             results.append(i**2)
    #         return results
    #     args = [x for x in range(100)]
    #     results = multiprocess.main(function=test_function, 
    #                                 args=args, 
    #                                 debug=True, jobid=None, cpus=cpus)
        
    #     self.assertEqual(len(results), 100)

    # def tearDown(self):
    #     [path.unlink() for path in self.testdir.iterdir() if 'bed1' not in path.name and 'bed2' not in path.name and 'test_' not in path.name]

class TestMain(unittest.TestCase):
    def setUp(self):
        self.srcdir = Path(__file__).parent.parent.parent
        self.testdir = self.srcdir / 'test_files'
        self.bed1 = [self.testdir / 'SRR1105736.tfit_bidirs.chr22.bed', 
                    self.testdir / 'SRR1105737.tfit_bidirs.chr22.bed']
        # self.bed1 = ' '.join([p.as_posix() for p in self.bed1])
        self.bed2 = [self.testdir / 'SRR1105738.tfit_bidirs.chr22.bed', 
                    self.testdir / 'SRR1105739.tfit_bidirs.chr22.bed']
        # self.bed2 = ' '.join([p.as_posix() for p in self.bed2])
        self.smallwindow = 150
        self.largewindow = 1500
        self.fimo_thresh = 0.000001
        self.bam1 = [self.testdir / 'SRR1105736.sorted.chr22.subsample.bam', 
                    self.testdir / 'SRR1105737.sorted.chr22.subsample.bam']
        # self.bam1 = ' '.join([p.as_posix() for p in self.bam1])
        self.bam2 = [self.testdir / 'SRR1105738.sorted.chr22.subsample.bam', 
                    self.testdir / 'SRR1105739.sorted.chr22.subsample.bam']
        # self.bam2 = ' '.join([p.as_posix() for p in self.bam2])
        self.label1 = 'condition1'
        self.label2 = 'condition2'
        self.genomefasta = self.testdir / 'chr22.fa'
        self.ranked_file = self.testdir / 'ranked_file.bed'

        self.count_file = self.testdir / 'count_file.header.bed'
        self.fimo_motifs = self.testdir / 'test_database.meme'
    
    def test_main(self):
        TFEA_path = Path(__file__).absolute().parent.parent
        command = ['python3', TFEA_path, '--output', self.testdir / 'test_output', 
                    '--bed1'] + self.bed1 + ['--bed2'] +  self.bed2 + ['--bam1'] + self.bam1 + ['--bam2'] + self.bam2 + ['--label1', self.label1, 
                    '--label2', self.label2,
                    '--genomefasta', self.genomefasta,
                    '--fimo_motifs', self.fimo_motifs]
        try:
            subprocess.check_output(command, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())
    
    def test_1rep(self):
        TFEA_path = Path(__file__).absolute().parent.parent
        command = ['python3', TFEA_path, 
                    '--output', self.testdir / 'test_output', 
                    '--bed1', self.bed1[0], 
                    '--bed2', self.bed2[0], 
                    '--bam1', self.bam1[0], 
                    '--bam2', self.bam2[0], 
                    '--label1', self.label1, 
                    '--label2', self.label2,
                    '--genomefasta', self.genomefasta,
                    '--fimo_motifs', self.fimo_motifs]
        try:
            subprocess.check_output(command, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())

    def test_html(self):
        TFEA_path = Path(__file__).absolute().parent.parent
        command = ['python3', TFEA_path, 
                    '--output', self.testdir / 'test_output', 
                    '--bed1', self.bed1[0], 
                    '--bed2', self.bed2[0], 
                    '--bam1', self.bam1[0], 
                    '--bam2', self.bam2[0], 
                    '--label1', self.label1, 
                    '--label2', self.label2,
                    '--genomefasta', self.genomefasta,
                    '--fimo_motifs', self.fimo_motifs, 
                    '--outputtype', 'html']
        try:
            subprocess.check_output(command, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise exceptions.SubprocessError(e.stderr.decode())


    # def tearDown(self):
    #     [path.unlink() for path in self.testdir.iterdir() if 'bed1' not in path.name and 'bed2' not in path.name and 'test_' not in path.name]

if __name__ == '__main__':
    unittest.main()