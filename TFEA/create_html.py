__author__ = 'Jonathan Rubin'

from config import *
import datetime

def createTFtext(TFresults,outputdir):
    TFresults = sorted(TFresults, key=lambda x: x[3])
    outfile = open(outputdir + 'results.txt', 'w')
    outfile.write('TF-Motif\tES\tNES\tP-value\tFDR\n')
    for val in TFresults:
        outfile.write('\t'.join([str(val[i]) for i in range(len(val))]) +  '\n')
    outfile.close()


def run(TFresults,outputdir,COMBINEtime,COUNTtime,DESEQtime,CALCULATEtime):
    #Creates results.txt which is a tab-delimited text file with the results    
    ##TFresults = sorted(TFresults, key=lambda x: x[3])
    outfile = open(outputdir + 'results.txt', 'w')
    outfile.write('TF-Motif\tES\tNES\tP-value\tFDR\n')
    for val in TFresults:
        outfile.write('\t'.join([str(val[i]) for i in range(len(val))]) +  '\n')
    outfile.close()

    #summary.html contains all user-defined variables, and also information about module used
    outfile = open(outputdir+'summary.html','w')
    outfile.write("""<!DOCTYPE html>
            <html>
            <head>
            <title>Variables Used</title>
            </head>
            <body>
                <h1>Variables Used</h1>
                <p>BEDS = """+str(BEDS)+"""</p>
                <p>LABEL1 = """+LABEL1+"""</p>
                <p>LABEL2 = """+LABEL2+"""</p>
                <p>BAM1 = """+str(BAM1)+"""</p>
                <p>BAM2 = """+str(BAM2)+"""</p>
                <p>SINGLEMOTIF = """+str(SINGLEMOTIF)+"""</p>
                <p>MOTIF_HITS = """+str(MOTIF_HITS)+"""</p>
                <p>OUTPUT = """+OUTPUT+"""
            </body>""")

    #For each TF motif with an FDR value less than a cutoff, an html file is created to be used in results.html
    for i in range(len(TFresults)):
        MOTIF_FILE,ES,NES,PVAL,POS,NEG,FDR = TFresults[i]
        positivelist = [x[0] for x in TFresults if x[2] > 0]
        negativelist = [x[0] for x in TFresults if x[2] < 0]
        
        if NES > 0:
            try:
                NEXT_MOTIF = positivelist[positivelist.index(MOTIF_FILE)+1]
            except IndexError:
                NEXT_MOTIF = positivelist[0]
            try:
                PREV_MOTIF = positivelist[positivelist.index(MOTIF_FILE)-1]
            except IndexError:
                PREV_MOTIF = positivelist[len(positivelist)]
        else:
            try:
                NEXT_MOTIF = negativelist[negativelist.index(MOTIF_FILE)+1]
            except IndexError:
                NEXT_MOTIF = negativelist[0]
            try:
                PREV_MOTIF = negativelist[negativelist.index(MOTIF_FILE)-1]
            except IndexError:
                PREV_MOTIF = negativelist[len(negativelist)]
        # if PVAL < FDRCUTOFF:
        outfile = open(outputdir + 'plots/' + MOTIF_FILE + '.results.html','w')
        outfile.write("""<!DOCTYPE html>
        <html>
        <head>
        <title>"""+MOTIF_FILE+""" Results</title>
        <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
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
        <body style="width: 1300px; overflow:scroll">
        <div style="float:left">
            <a href="./"""+PREV_MOTIF+""".results.html">PREV</a>
        </div>
        <div style="float:right">
            <a href="./"""+NEXT_MOTIF+""".results.html">NEXT</a>
        </div>
        <div style="text-align:center">
            <a href="../results.html">ALL</a>
        </div>
            <h1>"""+MOTIF_FILE+""" Results</h1>
        <div>
            <div style="float: middle; width: 1300px; overflow:scroll; padding-bottom:25px; padding-top:25px">
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>ES</th> 
                        <th>NES</th>
                        <th>P-value</th>
                        <th>FDR</th>
                        <th>Total Hits</th>
                        <th>Positive Hits</th> 
                        <th>Negative Hits</th>
                    </tr>
                    <tr>
                        <td>"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str("%.4g" % FDR)+"""</td>
                        <td>"""+str(int(POS+NEG))+"""</td>
                        <td>"""+str(int(POS))+"""</td>
                        <td>"""+str(int(NEG))+"""</td>
                    </tr>
                </table>
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; padding-top:50px">
                <img src="./"""+MOTIF_FILE+"""_enrichment_plot.png" alt="Enrichment Plot">
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; padding-top:50px">
                <img src="./"""+MOTIF_FILE+"""_meta_eRNA.png" alt="Meta Plot">
            </div>
        </div>
        <div>
            <div style="float: right; width: 600px; overflow: scroll">
                <p>Forward:</p>
                <img src="./"""+MOTIF_FILE.split('HO_')[1]+"""_direct.png" alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./"""+MOTIF_FILE.split('HO_')[1]+"""_revcomp.png" alt="Reverse Logo">
            </div>
            <div style="float:left; width: 600px overflow:scroll">
                <img src="./"""+MOTIF_FILE+"""_simulation_plot.png" alt="Simulation Plot">
            </div>
        </div>
        
        </body>
        </html>""")
        outfile.close()
        PREV_MOTIF = MOTIF_FILE

    outfile = open(outputdir+'results.html','w')
    outfile.write("""<!DOCTYPE html>
    <html>
    <head>
    <title>TFEA Results</title>
    <style>
    table {
        font-family: arial, sans-serif;
        border-collapse: collapse;
        width: 100%;
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
    <body style="width: 1300px; overflow:scroll">

    <h1>TFEA Results """ +LABEL1+ """ vs. """ +LABEL2+ """</h1>
    <div style="width: 1300px">
        <div style="float: left; width: 600px; padding-bottom: 50px">
            <img src="./plots/TFEA_NES_MA_Plot.png" alt="NES MA-Plot">
        </div>
        <div style="float: right; width:650px; padding-bottom: 50px">
            <img src="./plots/TFEA_Results_Moustache_Plot.png" alt="Moustache Plot (FDR vs. NES)">
        </div>
        <div id="Summary of Variables Used" style="float: right; width: 600px; padding-bottom: 315px">
            <p><a href="./Summary.html">Full Summary of Variables Used</a></p>
            <p><b>FDR < """ + str(FDRCUTOFF) + """</b></p>
            <table>
                <tr>
                    <th>Module</th>
                    <th>Switch</th>
                    <th>Time (hh:mm:ss)</th>
                </tr>
                <tr>
                    <td>COMBINE</td>
                    <td>"""+str(COMBINE)+"""</td>
                    <td>"""+str(datetime.timedelta(seconds=int(COMBINEtime)))+"""</td>
                </tr>
                <tr>
                    <td>COUNT</td>
                    <td>"""+str(COUNT)+"""</td>
                    <td>"""+str(datetime.timedelta(seconds=int(COUNTtime)))+"""</td>
                </tr>
                <tr>
                    <td>DESEQ</td>
                    <td>"""+str(DESEQ)+"""</td>
                    <td>"""+str(datetime.timedelta(seconds=int(DESEQtime)))+"""</td>
                </tr>
                <tr>
                    <td>CALCULATE</td>
                    <td>"""+str(CALCULATE)+"""</td>
                    <td>"""+str(datetime.timedelta(seconds=int(CALCULATEtime)))+"""</td>
                </tr>
            </table>   
        </div>
    </div>
    <div style="width: 1300px">
        <div style="float: left; width: 650px; overflow: scroll">
            <img src="./plots/TFEA_Pval_Histogram.png" alt="P-value Histogram">
        </div>
    </div>
    <div>
        <div id="Positive Enrichment Score" style="float: left; width: 600px; overflow:scroll">
            <h1>Positive Enrichment Score</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                    <th>Hits</th>
                </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,POS,NEG,FDR in TFresults:
        if NES > 0:
            if PVAL < FDR:
                outfile.write("""                <tr style="color: red;">
                        <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.3g" % PVAL)+"""</td>
                        <td>"""+str("%.3g" % FDR)+"""</td>
                        <td>"""+str(int(POS+NEG))+"""</td>
                    </tr>
                    """)
            else:
                outfile.write("""                <tr>
                        <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.3g" % PVAL)+"""</td>
                        <td>"""+str("%.3g" % FDR)+"""</td>
                        <td>"""+str(int(POS+NEG))+"""</td>
                    </tr>
                    """)


    outfile.write("""            </table>
        </div>

        <div id="Negative Enrichment Score" style="float: right; width: 600px; overflow:scroll">
            <h1>Negative Enrichment Score</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                    <th>Hits</th>
                </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,POS,NEG,FDR in TFresults:
        if NES < 0:
            if PVAL < FDR:
                outfile.write("""                <tr style="color: red;">
                        <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.3g" % PVAL)+"""</td>
                        <td>"""+str("%.3g" % FDR)+"""</td>
                        <td>"""+str(int(POS+NEG))+"""</td>
                    </tr>
                    """)
            else:
                outfile.write("""                <tr>
                        <td><a href="./plots/"""+MOTIF_FILE+""".results.html">"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.3g" % PVAL)+"""</td>
                        <td>"""+str("%.3g" % FDR)+"""</td>
                        <td>"""+str(int(POS+NEG))+"""</td>
                    </tr>
                    """)

    outfile.write("""        </table>
        </div>
    </div>

    </body>
    </html>""")

    outfile.close()

def single_motif(results,outputdir):
    MOTIF_FILE,ES,NES,PVAL,POS,NEG = results
    outfile = open(outputdir + MOTIF_FILE + '.results.html','w')
    outfile.write("""<!DOCTYPE html>
        <html>
        <head>
        <title>"""+MOTIF_FILE+""" Results</title>
        <style>
        table {
            font-family: arial, sans-serif;
            border-collapse: collapse;
            width: 100%;
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
        <body style="width: 1300px; overflow:scroll">
            <h1>"""+MOTIF_FILE+""" Results</h1>
        <div>
            <div style="float: middle; width: 1300px; overflow:scroll; padding-bottom:25px; padding-top:25px">
                <table> 
                    <tr>
                        <th>TF Motif</th>
                        <th>ES</th> 
                        <th>NES</th>
                        <th>P-value</th>
                        <th>Total Hits</th>
                        <th>Positive Hits</th> 
                        <th>Negative Hits</th>
                    </tr>
                    <tr>
                        <td>"""+MOTIF_FILE+"""</td>
                        <td>"""+str("%.3f" % ES)+"""</td>
                        <td>"""+str("%.3f" % NES)+"""</td>
                        <td>"""+str("%.4g" % PVAL)+"""</td>
                        <td>"""+str(int(POS+NEG))+"""</td>
                        <td>"""+str(int(POS))+"""</td>
                        <td>"""+str(int(NEG))+"""</td>
                    </tr>
                </table>
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; padding-top:50px">
                <img src="./plots/"""+MOTIF_FILE+"""_enrichment_plot.png" alt="Enrichment Plot">
            </div>
        </div>
        <div>
            <div style="float: left; width 1250px; padding-bottom:50px; padding-top:50px">
                <img src="./plots/"""+MOTIF_FILE+"""_meta_eRNA.png" alt="Meta Plot">
            </div>
        </div>
        <div>
            <div style="float: right; width: 600px; overflow: scroll">
                <p>Forward:</p>
                <img src="./plots/"""+MOTIF_FILE.split('HO_')[1]+"""_direct.png" alt="Forward Logo">
                <p></p>
                <p>Reverse:</p>
                <img src="./plots/"""+MOTIF_FILE.split('HO_')[1]+"""_revcomp.png" alt="Reverse Logo">
            </div>
            <div style="float:left; width: 600px overflow:scroll">
                <img src="./plots/"""+MOTIF_FILE+"""_simulation_plot.png" alt="Simulation Plot">
            </div>
        </div>
        
        </body>
        </html>""")
    outfile.close()


if __name__ == "__main__":
    TFresults = [['HO_P53_HUMAN.H10MO.B.bed',0.182143716966,6.22622338072,7.43595407471e-10,0.0]]
    filedir = './'
    figuredir = './'
    run(filedir,TFresults,figuredir)
