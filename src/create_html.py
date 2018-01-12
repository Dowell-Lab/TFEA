__author__ = 'Jonathan Rubin'

def run(filedir):
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
    <body style="width: 1100px; overflow:scroll">

    <h1>TFEA Results</h1>
    <div>
        <div id="Positively Enriched" style="float: left; width: 500px; overflow:scroll">
            <h1>Positively Enriched</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>
            </table>
        </div>

        <div id="Negatively Enriched" style="float: right; width: 500px; overflow:scroll">
            <h1>Negatively Enriched</h1>
            <table> 
                <tr>
                    <th>TF Motif</th>
                    <th>ES</th> 
                    <th>NES</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>
            </table>
        </div>
    </div>

    </body>
    </html>""")

if __name__ == "__main__":
    filedir = './'
    run(filedir)