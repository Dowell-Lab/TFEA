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
        self.bg1 = [self.testdir / 'SRR1105736.chr22.bedGraph', 
                    self.testdir / 'SRR1105737.chr22.bedGraph']
        self.bg2 = [self.testdir / 'SRR1105738.chr22.bedGraph', 
                    self.testdir / 'SRR1105739.chr22.bedGraph']
        self.label1 = 'DMSO'
        self.label2 = 'Nutlin'
        self.genomefasta = self.testdir / 'chr22.fa'
        self.ranked_file = self.testdir / 'ranked_file.bed'
        self.fasta_file = self.testdir / 'test_fasta_file.fa'

        self.count_file = self.testdir / 'count_file.header.bed'
        self.fimo_motifs = self.testdir / 'test_database.meme'
        # self.fimo_motifs = self.testdir.parent.parent.parent / 'motif_files' / 'best_curated_Human.meme'
        self.genomehits = self.testdir / 'test_genome_hits'

        # touch_command = "touch " + str(self.testdir / '*.bai')
        # subprocess.run(touch_command, check=True, shell=True)
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
                    '--metaprofile',
                    '--combine', 'mumerge',
                    # '--mdd_percent', '0.2',
                    '--md', '--mdd',
                    '--debug',
                    '--cpus', '2',
                    '--padjcutoff', '0.1']
        # try:
        print('\n============================================', 
                file=sys.stderr, flush=True)
        for output in execute(command):
            try:
                print('\t', output.decode(), end="", flush=True)
            except:
                print('\t', output, end="", flush=True)
            # subprocess.check_output(command, stderr=subprocess.PIPE)
        # except subprocess.CalledProcessError as e:
        #     raise exceptions.SubprocessError(e.stderr.decode())

        self.assertTrue((self.testdir / 'test_output' / 'results.txt').exists())
        self.assertTrue((self.testdir / 'test_output' / 'results.html').exists())
        # self.assertTrue((self.testdir / 'test_output' / 'summary.html').exists())
        # self.assertTrue((self.testdir / 'test_rep1' / 'md_results.txt').exists())
        # self.assertTrue((self.testdir / 'test_rep1' / 'mdd_results.txt').exists())

    # def test_fasta_file(self):
    #     shutil.rmtree(self.testdir / 'test_output', ignore_errors=True)
    #     TFEA_path = Path(__file__).absolute().parent.parent
    #     command = ['nice', '-n', '19',
    #                 'python3', TFEA_path, 
    #                 '--output', self.testdir / 'test_output', 
    #                 '--fasta_file', self.fasta_file, 
    #                 '--label1', self.label1, 
    #                 '--label2', self.label2,
    #                 '--fimo_motifs', self.fimo_motifs,  
    #                 '--motif_annotation', self.testdir / 'test_motif_annotation.bed',
    #                 '--output_type', 'html', 
    #                 '--plotall',
    #                 '--debug',
    #                 '--cpus', '2',
    #                 '--padjcutoff', '0.1']
    #     # try:
    #     print('\n============================================', 
    #             file=sys.stderr, flush=True)
    #     for output in execute(command):
    #         try:
    #             print('\t', output.decode(), end="", flush=True)
    #         except:
    #             print('\t', output, end="", flush=True)
    #         # subprocess.check_output(command, stderr=subprocess.PIPE)
    #     # except subprocess.CalledProcessError as e:
    #     #     raise exceptions.SubprocessError(e.stderr.decode())

    #     self.assertTrue((self.testdir / 'test_output' / 'results.txt').exists())
    #     self.assertTrue((self.testdir / 'test_output' / 'results.html').exists())
    #     # self.assertTrue((self.testdir / 'test_output' / 'summary.html').exists())
    #     # self.assertTrue((self.testdir / 'test_rep1' / 'md_results.txt').exists())
    #     # self.assertTrue((self.testdir / 'test_rep1' / 'mdd_results.txt').exists())

    # def test_genome_hits(self):
    #     shutil.rmtree(self.testdir / 'test_output', ignore_errors=True)
    #     TFEA_path = Path(__file__).absolute().parent.parent
    #     command = ['nice', '-n', '19',
    #                 'python3', TFEA_path, 
    #                 '--output', self.testdir / 'test_output', 
    #                 '--bed1', self.bed1[0], 
    #                 '--bed2', self.bed2[0], 
    #                 '--bam1', self.bam1[0], 
    #                 '--bam2', self.bam2[0], 
    #                 '--label1', self.label1, 
    #                 '--label2', self.label2,
    #                 '--genomefasta', self.genomefasta,
    #                 '--fimo_motifs', self.fimo_motifs,
    #                 '--scanner', 'genome hits',
    #                 '--genomehits', self.genomehits,  
    #                 '--motif_annotation', self.testdir / 'test_motif_annotation.bed',
    #                 '--output_type', 'html', 
    #                 '--plotall',
    #                 '--combine', 'merge all',
    #                 '--debug']
    #     # try:
    #     print('\n============================================')
    #     for output in execute(command):
    #         print('\t', output.decode(), end="", flush=True)
    #         # subprocess.check_output(command, stderr=subprocess.PIPE)
    #     # except subprocess.CalledProcessError as e:
    #     #     raise exceptions.SubprocessError(e.stderr.decode())

    #     self.assertTrue((self.testdir / 'test_output' / 'results.txt').exists())
    #     self.assertTrue((self.testdir / 'test_output' / 'results.html').exists())
    #     # self.assertTrue((self.testdir / 'test_output' / 'summary.html').exists())
    #     # self.assertTrue((self.testdir / 'test_rep1' / 'md_results.txt').exists())
    #     # self.assertTrue((self.testdir / 'test_rep1' / 'mdd_results.txt').exists())

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
    #         subprocess.check_output(command, stderr=subprocess.PIPE)
    #     # except subprocess.CalledProcessError as e:
    #     #     raise exceptions.SubprocessError(e.stderr.decode())

    #     self.assertTrue((self.testdir / 'test_rep2' / 'results.txt').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'results.html').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'summary.html').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'md_results.txt').exists())
    #     self.assertTrue((self.testdir / 'test_rep2' / 'mdd_results.txt').exists())

    # def test_combine(self):
    #     from TFEA import combine
    #     merged_bed1 = combine.merge_bed(beds=self.bed1+self.bed2)
    #     merged_bed2 = combine.merge_bed(beds=self.bed2+self.bed1)
        
    #     print(len(merged_bed1), len(merged_bed2))

if __name__ == '__main__':
    import sys
    from pathlib import Path
    # Add TFEA srcdirectory into path
    srcdirectory = Path(__file__).absolute().parent.parent
    sys.path.insert(0, srcdirectory)
    unittest.main(verbosity=2, failfast=True)