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
import math
import datetime
import numpy as np
import matplotlib.pyplot as plt

#Functions
#==============================================================================
def txt_output(results=None, outputdir=None, outname='results.txt', 
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
def html_output(results=None, module_list=None, columns=None):
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
    #Using a config file
    import config
    outputdir = config.OUTPUTDIR
    label1 = config.LABEL1
    label2 = config.LABEL2
    padj_cutoff = config.PADJCUTOFF
    plot = config.PLOTALL

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
                <img src="./plots/TFEA_NES_MA_Plot.png" alt="NES MA-Plot">
            </div>
            <div style="float: left; width: 45%">
                <img src="./plots/DESEQ_MA_Plot.png" alt="DE-Seq MA-Plot">
            </div>
        </div>
        <div class="row">
            <div style="float: left; width:45%">
                <img src="./plots/TFEA_Pval_Histogram.png" alt="TFEA P-Value \
                    Histogram">
            </div>
            <div id="Summary of Variables Used" style="float: right; \
                width: 45%">
                <p><a href="./Summary.html">Full Summary of Variables Used\
                </a></p>
                <p><b>PADJ < """ + str(padj_cutoff) + """</b></p>
                <table>
                    <tr>
                        <th>Module</th>
                        <th>Value</th>
                        <th>Time (hh:mm:ss)</th>
                    </tr>""")
    for (module, value, time) in module_list:
        outfile.write("""<tr>
                        <td>"""+module+"""</td>
                        <td>"""+str(value)+"""</td>
                        <td>"""+str(datetime.timedelta(
                                    seconds=int(time)))
                        +"""</td>
                    </tr>""")
    outfile.write("""
                </table>   
            </div>
        </div>
        <div>
            <div id="Positive Enrichment Score" style="float: left; width:45%">
                <h1>Positive Enrichment Score</h1>
                <table> 
                    <tr>
    """)

    for column_name in columns:
        outfile.write("""                        <th>"""+column_name+"""</th>
    """)

    outfile.write("""</tr>
                """)

    for motif, enrichment, hits, p_val, p_adj in results:
        if enrichment > 0:
            if p_adj < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p_val)+"""</td>
                <td>"""+str("%.3g" % p_adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p_val)+"""</td>
                <td>"""+str("%.3g" % p_adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)

            else:
                outfile.write("""
            <tr>
                <td>"""+motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p_val)+"""</td>
                <td>"""+str("%.3g" % p_adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)


    outfile.write("""            
        </table>
    </div>

    <div id="Negative Enrichment Score" style="float: right; width: 45%">
        <h1>Negative Enrichment Score</h1>
        <table> 
            <tr>
                <th>TF Motif</th>
                <th>ENRICHMENT("""+config.ENRICHMENT+""")</th>
                <th>P-value</th>
                <th>PADJ</th>
                <th>HITS</th>
            </tr>
                """)

    for [motif, enrichment, hits, p_val, p_adj] in results:
        if enrichment < 0:
            if p_adj < padj_cutoff:
                outfile.write("""
            <tr style="color: red;">
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p_val)+"""</td>
                <td>"""+str("%.3g" % p_adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)
            elif plot:
                outfile.write("""
            <tr>
                <td><a href="./plots/"""+motif+""".results.html">"""
                    +motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p_val)+"""</td>
                <td>"""+str("%.3g" % p_adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
                    """)
            else:
                outfile.write("""
            <tr>
                <td>"""+motif+"""</td>
                <td>"""+str("%.4f" % enrichment)+"""</td>
                <td>"""+str("%.3g" % p_val)+"""</td>
                <td>"""+str("%.3g" % p_adj)+"""</td>
                <td>"""+str(hits)+"""</td>
            </tr>
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
    with open(os.path.join(outputdir,'summary.html'),'w') as outfile:
        outfile.write("""<!DOCTYPE html>
                <html>
                <head>
                <title>Variables Used</title>
                </head>
                <body>
                    <h1>Variables Used</h1>""")
        for key in config_object:
            for item in config_object[key]:
                outfile.write('<p>' + item.upper() + '=' 
                                + config_object[key][item] + '</p>\n')
        outfile.write('</body>')

#==============================================================================
def create_motif_result_html(results=None):
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
    #Using a config file
    import config
    outputdir = config.OUTPUTDIR
    padj_cutoff = config.PADJCUTOFF
    singlemotif = config.SINGLEMOTIF
    plot = config.PLOTALL

    #For each TF motif with an PADJ value less than a cutoff, an html file is 
    #created to be used in results.html
    positivelist = [x[0] for x in results 
                    if x[2] > 0 and (plot or x[-1] < padj_cutoff)]
    negativelist = [x[0] for x in results 
                    if x[2] < 0 and (plot or x[-1] < padj_cutoff)]

    for i in range(len(results)):
        motif, enrichment, hits, p_val, p_adj = results[i] 
        if plot or p_adj < padj_cutoff or singlemotif:
            if enrichment > 0:
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
            direct_logo = motif.strip('HO_') + "_direct.png"
            reverse_logo = motif.strip('HO_') + "_revcomp.png"
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
                        <th>TF Motif</th>
                        <th>AUC</th> 
                        <th>P-value</th>
                        <th>PADJ</th>
                        <th>HITS</th>
                    </tr>
                    <tr>
                        <td>"""+motif+"""</td>
                        <td>"""+str("%.2f" % enrichment)+"""</td>
                        <td>"""+str("%.4g" % p_val)+"""</td>
                        <td>"""+str("%.4g" % p_adj)+"""</td>
                        <td>"""+str(hits)+"""</td>
                    </tr>
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
            <div style="float: right; width: 600px">
                <p>Forward:</p>
                <img src="./"""+direct_logo+"""" \
                    alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./"""+reverse_logo+"""" \
                    alt="Reverse Logo">
            </div>
            <div style="float:left; width: 600px">
                <img src="./"""+motif+"""_simulation_plot.png" \
                    alt="Simulation Plot">
            </div>
        </div>

    </body>
    </html>""")
            outfile.close()
            PREV_MOTIF = motif

#==============================================================================
def plot_deseq_MA(deseq_file=None, label1=None, label2=None, figuredir=None):
    '''Plots the DE-Seq MA-plot using the full regions of interest and saves it
    to the figuredir directory created in TFEA output folder

    Parameters
    ----------
    deseqfile : string
        full path to the deseq file (specifically .res.txt)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    figuredir : string
        full path to figure directory in output directory (created by TFEA)

    Returns
    -------
    None
    '''
    up_x = list()
    up_y = list()
    up_p = list()
    dn_x = list()
    dn_y = list()
    dn_p = list()
    with open(deseq_file,'r') as F:
        header = F.readline().strip('\n').split('\t')
        basemean_index = header.index('"baseMean"')
        log2fc_index = header.index('"log2FoldChange"')
        for line in F:
            line = line.strip('\n').split('\t')
            try:
                log2fc = float(line[log2fc_index+1])
                basemean = math.log(float(line[basemean_index+1]),10)
                pval = float(line[-2])
                if log2fc > 0:
                    up_x.append(basemean)
                    up_y.append(log2fc)
                    up_p.append(pval)
                else:
                    dn_x.append(basemean)
                    dn_y.append(log2fc)
                    dn_p.append(pval)
            except:
                pass

    x = [x for _,x in sorted(zip(up_p,up_x))] \
        + [x for _,x in sorted(zip(dn_p,dn_x),reverse=True)]

    y = [y for _,y in sorted(zip(up_p,up_y))] \
        + [y for _,y in sorted(zip(dn_p,dn_y),reverse=True)]

    c = plt.cm.RdYlGn(np.linspace(0, 1, len(x)))

    #Creates an MA-Plot of the region expression
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    plt.scatter(x=x,y=y,color=c,edgecolor='')
    ax.set_title("DE-Seq MA-Plot",fontsize=14)
    ax.set_ylabel("Log2 Fold-Change ("+label2+"/"+label1+")",fontsize=14)
    ax.set_xlabel("Log10 Average Expression",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', 
                    labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', 
                    labelbottom='on')
    plt.savefig(os.path.join(figuredir, 'DESEQ_MA_Plot.png'),
                bbox_inches='tight')