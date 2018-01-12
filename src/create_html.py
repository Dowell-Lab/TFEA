__author__ = 'Jonathan Rubin'

def run(filedir,TFresults,figuredir):
    outfile = open(filedir+'results.html','w')
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
    <body style="width: 1250px; overflow:scroll">

    <h1>TFEA Results</h1>
    <div>
        <div id="Positively Enriched" style="float: left; width: 600px; overflow:scroll">
            <h1>Positively Enriched</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,FDR in TFresults:
        if NES > 0:
            outfile.write("""                <tr>
                """)
            outfile.write("""                    <td>"""+MOTIF_FILE+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.3f" % ES)+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.3f" % NES)+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.4g" % PVAL)+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.4g" % FDR)+"""</td>
                """)
            outfile.write("""                </tr>
                """)


    outfile.write("""            </table>
        </div>

        <div id="Negatively Enriched" style="float: right; width: 600px; overflow:scroll">
            <h1>Negatively Enriched</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>
                """)

    for MOTIF_FILE,ES,NES,PVAL,FDR in TFresults:
        if NES < 0:
            outfile.write("""                <tr>
                """)
            outfile.write("""                    <td>"""+MOTIF_FILE+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.3f" % ES)+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.3f" % NES)+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.4g" % PVAL)+"""</td>
                """)
            outfile.write("""                    <td>"""+str("%.4g" % FDR)+"""</td>
                """)
            outfile.write("""                </tr>
                """)

    outfile.write("""        </table>
        </div>
    </div>

    </body>
    </html>""")

if __name__ == "__main__":
    TFresults = [['HO_P53_HUMAN.H10MO.B.bed',0.182143716966,6.22622338072,7.43595407471e-10,0.0]]
    filedir = './'
    figuredir = './'
    run(filedir,TFresults,figuredir)