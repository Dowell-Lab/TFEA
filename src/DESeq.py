__author__ = 'Jonathan Rubin'

import os

def write_script(LABEL1,LABEL2,BAM1,BAM2,scriptdir,filedir,count_file):
    outfile = open(scriptdir + 'DESeq.R','w')
    outfile.write('sink("'+filedir+'DESeq.Rout")\n')
    outfile.write('library("DESeq")\n')
    outfile.write('data <- read.delim("'+count_file+'", sep="\t", header=TRUE)\n')
    outfile.write('countsTable <- subset(data, select=c('+', '.join([str(i) for i in range(5,5+len(BAM1)+len(BAM2))])+'))\n')
    outfile.write('rownames(countsTable) <- data$region\n')
    outfile.write('conds <- c(' + ', '.join(['"'+LABEL1+'"']*len(BAM1)) + ', ' + ', '.join(['"'+LABEL2+'"']*len(BAM2)) + ')\n')
    outfile.write('cds <- newCountDataSet( countsTable, conds )\n')
    outfile.write('cds <- estimateSizeFactors( cds )\n')
    outfile.write('sizeFactors(cds)\n')
    outfile.write('cds <- estimateDispersions( cds )\n')
    outfile.write('res <- nbinomTest( cds, "'+LABEL1+'", "'+LABEL2+'" )\n')
    outfile.write('write.table(res, file = "'+filedir+'DESeq.res.txt", append = FALSE, sep= "\t" )\n')
    outfile.write('sink()')

def run(LABEL1,LABEL2,BAM1,BAM2,scriptdir,filedir):
    write_script(LABEL1,LABEL2,BAM1,BAM2,scriptdir,filedir)
    os.system("R < " + scriptdir + "DESeq.R --no-save")
    return filedir+'DESeq.res.txt'
