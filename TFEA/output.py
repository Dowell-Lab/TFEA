#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This file contains a list of functions associated with outputting a ranked
    bed file based on some metric.
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
import math
import time
import datetime
import numpy as np

from TFEA import multiprocess
from TFEA import plot

#Main Script
#==============================================================================
def main(use_config=True, outputdir=None, results=None, md_results=None, 
            mdd_results=None, motif_distances=None, md=False, mdd=False, 
            debug=False, jobid=None, output_type=False, p_cutoff=None, 
            dpi=100):
    '''This script creates output files associated with TFEA
    '''
    start_time = time.time()
    if use_config:
        from TFEA import config
        outputdir = config.vars['OUTPUT']
        figuredir = config.vars['FIGUREDIR']
        results=config.vars['RESULTS']
        md_results=config.vars['MD_RESULTS']
        mdd_results=config.vars['MDD_RESULTS']
        motif_distances = config.vars['MOTIF_DISTANCES']
        md = config.vars['MD']
        mdd = config.vars['MDD']
        debug = config.vars['DEBUG']
        jobid = config.vars['JOBID']
        output_type = config.vars['OUTPUT_TYPE']
        p_cutoff = config.vars['PADJCUTOFF']
        padj_cutoff = config.vars['PADJCUTOFF']
        label1 = config.vars['LABEL1']
        label2 = config.vars['LABEL2']
        plotall = config.vars['PLOTALL']
        singlemotif = config.vars['SINGLEMOTIF']
        dpi = config.vars['DPI']


    print("Creating output...", end=' ', flush=True, file=sys.stderr)
    TFEA_header = ['#TF', 'AUC', 'Events', 'Z-score','STD','FPKM', 'P-adj']
    sort_index = [3, 2, 1, -1]
    txt_output(outputdir=outputdir, results=results, outname='results.txt', 
                sortindex=sort_index, header=TFEA_header)
    plot.plot_global_MA(results, p_cutoff=p_cutoff, title='TFEA MA-Plot', 
                    xlabel='Log10(Motif Hits)', 
                    ylabel='Area Under the Curve (AUC)', 
                    savepath=figuredir / 'TFEA_MA.png', 
                    dpi=dpi)
    plot.plot_global_volcano(results, p_cutoff=p_cutoff, title='TFEA Volcano Plot', 
                    xlabel='Area Under the Curve (AUC)', 
                    ylabel='-log10(P-adj)', 
                    savepath=figuredir / 'TFEA_volcano.png', 
                    dpi=dpi)
    if md:
        header = ['#TF', 'MD-Score', 'Events', 'p-val']
        txt_output(outputdir=outputdir, results=md_results, 
                    outname='md_results.txt', header=header, sortindex=[-1])
        plot.plot_global_MA(md_results, p_cutoff=p_cutoff, 
                                title='MD MA-Plot', 
                                xlabel='Log10(Motif Hits)', 
                                ylabel='MD-Score Difference', 
                                savepath=figuredir / 'MD_MA.png',  
                                dpi=dpi)
        plot.plot_global_volcano(md_results, p_cutoff=p_cutoff, 
                                    title='MD Volcano Plot', 
                                    xlabel='MD-Score Difference', 
                                    ylabel='-log10(P-val)', 
                                    savepath=figuredir / 'MD_volcano.png', 
                                    dpi=dpi)
    if mdd:
        header = ['#TF', 'MD-Score', 'Events', 'p-val']
        txt_output(outputdir=outputdir, results=mdd_results, 
                    outname='mdd_results.txt', header=header, sortindex=[-1])
        plot.plot_global_MA(mdd_results, p_cutoff=p_cutoff, 
                                title='MDD MA-Plot', 
                                xlabel='Log10(Motif Hits)', 
                                ylabel='Differential MD-Score Difference', 
                                savepath=figuredir / 'MDD_MA.png', 
                                dpi=dpi)
        plot.plot_global_volcano(mdd_results, p_cutoff=p_cutoff, 
                                    title='MDD Volcano Plot', 
                                    xlabel='Differential MD-Score Difference', 
                                    ylabel='-log10(P-val)', 
                                    savepath=figuredir / 'MDD_volcano.png', 
                                    dpi=dpi)
    
    total_time = time.time() - start_time
    if use_config:
        config.vars['OUTPUTtime'] = total_time
    if output_type == 'html':
        if use_config:
            # summary_html_output(config_object=config.vars, outputdir=outputdir)
            module_list = [('COMBINE', config.vars['COMBINE'], config.vars['COMBINEtime']),
                        ('RANK', config.vars['RANK'], config.vars['RANKtime']), 
                        ('SCANNER', config.vars['SCANNER'], config.vars['SCANNERtime']), 
                        ('ENRICHMENT', config.vars['ENRICHMENT'], config.vars['ENRICHMENTtime']), 
                        ('OUTPUT', config.vars['OUTPUT_TYPE'], config.vars['OUTPUTtime'])]
        else:
            module_list = []
        create_motif_result_htmls(results=results, results_header=TFEA_header, 
                                    outputdir=outputdir, 
                                    padj_cutoff=padj_cutoff, 
                                    singlemotif=singlemotif, 
                                    plotall=plotall, auc_index=1, 
                                    padj_index=-1)
        html_output(results=results, results_header=TFEA_header,
                    module_list=module_list, 
                    outputdir=outputdir, label1=label1, label2=label2, 
                    padj_cutoff=padj_cutoff, plotall=plotall, auc_index=1, 
                    padj_index=-1, sortindex=sort_index)
        
    print("done in: " + str(datetime.timedelta(seconds=int(total_time))), file=sys.stderr)

    if debug:
        multiprocess.current_mem_usage(jobid)


#Functions
#==============================================================================
def txt_output(results=None, outputdir=None, outname=None, 
                header=None, sortindex=None):
    with open(os.path.join(outputdir, outname), 'w') as outfile:
        if type(header) == list:
            outfile.write('\t'.join(header) + '\n')
        elif type(header) == str:
            outfile.write(header + '\n')
        for index in sortindex:
            results.sort(key=lambda x: x[index])
        for values in results:
            outfile.write('\t'.join([str(x) for x in values])+'\n')

#==============================================================================
def html_output(results=None, module_list=None, outputdir=None,
                label1=None, label2=None, padj_cutoff=None, plotall=None, 
                results_header=None, auc_index=1, padj_index=-1, sortindex=None):
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plot : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    for index in sortindex:
        results.sort(key=lambda x: x[index])
    outfile = open(os.path.join(outputdir, 'results.html'),'w')
    outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>TFEA Results """ + label1 + """ vs. """ + label2 +"""</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">

        <h1>TFEA Results """ + label1 + """ vs. """ + label2 + """</h1>
        <div class="row">
            <div style="float: left; width: 45%">
                <img src="./plots/TFEA_MA.png" alt="TFEA MA-Plot">
            </div>
            <div style="float: left; width: 45%">
                <img src="./plots/TFEA_volcano.png" alt="DE-Seq MA-Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: left; width:45%">
                <img src="./plots/DESEQ_MA_Plot.png" alt="Differential MD-Score P-Value \
                    Histogram">
            </div>
            <div id="User Inputs" style="float: right; width: 45%">
                <p><b>PADJ < """ + str(padj_cutoff) + """</b></p>
                <p><a href="./inputs.txt">User Inputs</a></p>
                <p><a href="./plots/MD_MA.png">MD MA-Plot</a> | \
                    <a href="./plots/MD_volcano.png">MD VolcanoPlot</a></p>
                <p><a href="./plots/MDD_MA.png">MDD MA-Plot</a> | \
                    <a href="./plots/MDD_volcano.png">MDD VolcanoPlot</a></p>
                <table>
                    <tr>
                        <th>Module</th>
                        <th>Value</th>
                        <th>Time (hh:mm:ss)</th>
                    </tr>
                """)
    total_time = 0
    if len(module_list) != 0:
        for module, value, time in module_list:
            total_time += time
            outfile.write("""<tr>
                        <td>"""+module+"""</td>
                        <td>"""+str(value)+"""</td>
                        <td>"""+str(datetime.timedelta(seconds=int(time)))
                        +"""</td>
                    </tr>""")
        outfile.write("""<tr>
                        <td><b>Total</b></td>
                        <td> </td>
                        <td><b>"""+str(datetime.timedelta(seconds=int(total_time)))
                        +"""</b></td>
                    </tr>""")
    outfile.write("""
                </table>   
            </div>
        </div>
        <div>
            <div id="Positive AUC Value" style="float: left; width:45%">
                <h1>Positive AUC Value</h1>
                <table> 
                    <tr>
                    """)
    for label in results_header:
        outfile.write("<th>" + label + "</th>\n")

    for motif_result in results:
        auc = motif_result[auc_index]
        p_adj = motif_result[padj_index]
        motif = motif_result[0]
        if auc >= 0:
            if p_adj < padj_cutoff:
                motif = motif_result[0]
                outfile.write("""
                    </tr>
            <tr style="color: red;">
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>""")
                for number_result in motif_result[1:]:
                    try:
                        outfile.write("<td>" + str("%.3g" % number_result) + "</td>\n")
                    except TypeError:
                        outfile.write("<td>" + str(number_result) + "</td>\n")
                outfile.write("""            </tr>
                    """)
            elif plotall:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>""")
                for number_result in motif_result[1:]:
                    try:
                        outfile.write("<td>" + str("%.3g" % number_result) + "</td>\n")
                    except TypeError:
                        outfile.write("<td>" + str(number_result) + "</td>\n")
                outfile.write("""            </tr>
                    """)

            else:
                outfile.write("""
            <tr>
                <td>"""+motif+"""</td>""")
                for number_result in motif_result[1:]:
                    try:
                        outfile.write("<td>" + str("%.3g" % number_result) + "</td>\n")
                    except TypeError:
                        outfile.write("<td>" + str(number_result) + "</td>\n")
                outfile.write("""            </tr>
                    """)


    outfile.write("""            
        </table>
    </div>

    <div id="Negative AUC Value" style="float: right; width: 45%">
        <h1>Negative AUC Value</h1>
        <table> 
            <tr>
                """)
    for label in results_header:
        outfile.write("<th>" + label + "</th>\n")

    for motif_result in results:
        auc = motif_result[auc_index]
        p_adj = motif_result[padj_index]
        motif = motif_result[0]
        if auc < 0:
            if p_adj < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>""")
                for number_result in motif_result[1:]:
                    try:
                        outfile.write("<td>" + str("%.3g" % number_result) + "</td>\n")
                    except TypeError:
                        outfile.write("<td>" + str(number_result) + "</td>\n")
                outfile.write("""            </tr>
                    """)
            elif plotall:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>""")
                for number_result in motif_result[1:]:
                    try:
                        outfile.write("<td>" + str("%.3g" % number_result) + "</td>\n")
                    except TypeError:
                        outfile.write("<td>" + str(number_result) + "</td>\n")
                outfile.write("""            </tr>
                    """)
            else:
                outfile.write("""
            <tr>
                <td>"""+motif+"""</td>""")
                for number_result in motif_result[1:]:
                    try:
                        outfile.write("<td>" + str("%.3g" % number_result) + "</td>\n")
                    except TypeError:
                        outfile.write("<td>" + str(number_result) + "</td>\n")
                outfile.write("""            </tr>
                    """)

    outfile.write("""        
            </table>
        </div>
        </div>

    </body>
    </html>""")

    outfile.close()

#==============================================================================
def summary_html_output(config_object=None, outputdir=None):
    exclude = ['MOTIF_DISTANCES','MD_DISTANCES1', 'MD_DISTANCES2', 
                'MDD_DISTANCES1', 'MDD_DISTANCES2', 'PVALS', 'FCS', 
                'META_PROFILE', 'RESULTS', 'MD_RESULTS', 'MDD_RESULTS']
    with open(os.path.join(outputdir,'summary.html'),'w') as outfile:
        outfile.write("""<!DOCTYPE html>
                <html>
                <head>
                <title>Variables Used</title>
                </head>
                <body>
                    <a href="""+os.path.join(outputdir,'results.html')+""">BACK</a>
                    <h1>Variables Used</h1>""")
        for key in config_object:
            if key not in exclude:
                value = config_object[key]
                value = str(value)
                outfile.write('<p>' + key + '=' + value + '</p>\n')
        outfile.write('</body>')

#==============================================================================
def create_motif_result_htmls(results=None, outputdir=None, padj_cutoff=None,
                                singlemotif=None, plotall=None, 
                                results_header=None, auc_index=1, padj_index=3):
    '''Creates the main html output and also individual html outputs for each
        motif
    
    Parameters
    ----------
    outputdir : string
        the full path to the output directory created by TFEA

    beds : list or array
        a list of full paths to bed files to be considered as regions of 
        interest

    label1 : string
        an informative label describing sample corresponding to condition1

    label2 : string
        an informative label describing sample corresponding to condition2

    bam1 : list or array
        a list of full paths to bam files corresponding to condition1

    bam2 : list or array
        a list of full paths to bam files corresponding to condition2

    singlemotif : boolean or string
        either False if all motifs should be considered in TFEA or the name of
        a specific motif to be analyzed

    motif_hits : string
        the full path to a directory containing motif hits across the genome

    output : string
        the full path to a user-specified output directory. TFEA will create
        a new folder within this directory - this is called outputdir

    padj_cutoff : float
        the cutoff value for determining significance

    plotall : boolean
        a switch that controls whether all motifs are plotted or just 
        significant ones defined by the p-adj cutoff

    combine : boolean
        a switch that determines whether bed files within the beds variable
        get combined and merged using bedtools

    count : boolean
        a switch that controls whether reads are counted over the regions of
        interest

    deseq : boolean
        a switch that controls whether DE-Seq is performed on the inputted
        regions that have been counted over

    calculate : boolean
        a switch that determines whether the TFEA calculation is performed

    TFresults : list or array
        a list of lists contining 'enrichment' scores, normalized 'enrichment' 
        scores, p-value, p-adj, and number of hits for each individual motif

    COMBINEtime : float
        the time it took to combine and merge the bed files using bedtools

    COUNTtime : float
        the time it took to count reads over regions of interest

    DESEQtime : float
        the time it took to perform DE-Seq using the counts file

    CALCULATEtime : float
        the time it took to perform TFEA

    Returns
    -------
    None
    '''
    #For each TF motif with an PADJ value less than a cutoff, an html file is 
    #created to be used in results.html
    positivelist = [x[0] for x in results 
                    if x[1] >= 0 and (plotall or x[-1] < padj_cutoff)]
    negativelist = [x[0] for x in results 
                    if x[1] < 0 and (plotall or x[-1] < padj_cutoff)]

    for i in range(len(results)):
        motif = results[i][0]
        auc = results[i][auc_index]
        p_adj = results[i][padj_index]
        if plotall or p_adj < padj_cutoff: # or motif in singlemotif:
            if auc >= 0:
                try:
                    NEXT_MOTIF = positivelist[positivelist.index(motif)+1]
                except IndexError:
                    NEXT_MOTIF = positivelist[0]
                try:
                    PREV_MOTIF = positivelist[positivelist.index(motif)-1]
                except IndexError:
                    PREV_MOTIF = positivelist[len(positivelist)]
            else:
                try:
                    NEXT_MOTIF = negativelist[negativelist.index(motif)+1]
                except IndexError:
                    NEXT_MOTIF = negativelist[0]
                try:
                    PREV_MOTIF = negativelist[negativelist.index(motif)-1]
                except IndexError:
                    PREV_MOTIF = negativelist[len(negativelist)]
            direct_logo = "logo" + motif.replace('.','_') + ".png"
            reverse_logo = "logo_rc" + motif.replace('.','_') + ".png"
            outfile = open(os.path.join(outputdir, 'plots', motif 
                            + '.results.html'),'w')
            outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>"""+motif+""" Results</title>
    <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }

        .row {
        display: flex; /* equal height of the children */
        width: 100%;
        padding-bottom: 50px
        }

        img {
            max-width: 100%;
            max-height: 100%;
        }

        td, th {
            border: 1px solid #dddddd;
            text-align: left;
            padding: 8px;
        }

        tr:nth-child(even) {
            background-color: #dddddd;
        }
    </style>
    </head>
    <body style="width:1300px; margin:0 auto;">
        <div>
            <div style="float:left">
                <a href="./"""+PREV_MOTIF+""".results.html">PREV</a>
            </div>
            <div style="float:right">
                <a href="./"""+NEXT_MOTIF+""".results.html">NEXT</a>
            </div>
            <div style="text-align:center">
                <a href="../results.html">ALL</a>
        </div>
        <div class="row">
        </div>
            <h1>"""+motif+""" Results</h1>
        <div>
            <div style="float: middle; width: 1300px; padding-bottom:25px; \
                padding-top:25px">
                <table> 
                    <tr>
                """)
            for label in results_header:
                outfile.write("<th>" + label + "</th>\n")
            outfile.write("""
                    </tr>
                    <tr>
                        <td>"""+motif+"""</td>""")
            for number_result in results[i][1:]:
                outfile.write("<td>" + str("%.5g" % number_result) + "</td>\n")
            outfile.write("""            </tr>
                    
                </table>
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; \
                padding-top:50px">
                <img src="./"""+motif+"""_enrichment_plot.png" \
                    alt="Enrichment Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: left; width: 600px; padding-right:50px">
                <p>Forward:</p>
                <img src="./"""+direct_logo+"""" \
                    alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./"""+reverse_logo+"""" \
                    alt="Reverse Logo">
            </div>
            <div style="float:right; width: 600px">
                <img src="./"""+motif+"""_simulation_plot.png" \
                    alt="Simulation Plot">
            </div>
        </div>

    </body>
    </html>""")
            outfile.close()
            PREV_MOTIF = motif