__author__ = 'Jonathan Rubin'

import os
from config import LABEL1, LABEL2, BAM1, BAM2, FILEDIR, COUNT_FILE

def write_script():
    if (len(BAM1) > 1 and len(BAM2) > 1):
        outfile = open(FILEDIR + 'DESeq.R','w')
        outfile.write('sink("'+FILEDIR+'DESeq.Rout")\n')
        outfile.write('library("DESeq2")\n')
        outfile.write('data <- read.delim("'+COUNT_FILE+'", sep="\t", header=TRUE)\n')
        outfile.write('countsTable <- subset(data, select=c('+', '.join([str(i) for i in range(5,5+len(BAM1)+len(BAM2))])+'))\n')
        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- as.data.frame(c(' + ', '.join(['"'+LABEL1+'"']*len(BAM1)) + ', ' + ', '.join(['"'+LABEL2+'"']*len(BAM2)) + '))\n')
        outfile.write('colnames(conds) <- c("treatment")\n')
        outfile.write('ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsTable, colData = conds, design = ~ treatment)\n')
        outfile.write('dds <- DESeq(ddsFullCountTable)\n')
        outfile.write('res1 <- results(dds,alpha = 0.05, contrast=c("treatment","'+LABEL2+'","'+LABEL1+'"))\n')
        outfile.write('resShrink <- lfcShrink(dds, res = res1, contrast = c("treatment","'+LABEL2+'","'+LABEL1+'"))\n')
        outfile.write('resShrink$fc <- 2^(resShrink$log2FoldChange)\n')
        outfile.write('res <- resShrink[c(1:3,7,4:6)]\n')
        outfile.write('write.table(res, file = "'+FILEDIR+'DESeq.res.txt", append = FALSE, sep= "\t" )\n')
        outfile.write('sink()')
    else:
        outfile = open(FILEDIR + 'DESeq.R','w')
        outfile.write('sink("'+FILEDIR+'DESeq.Rout")\n')
        outfile.write('library("DESeq")\n')
        outfile.write('data <- read.delim("'+COUNT_FILE+'", sep="\t", header=TRUE)\n')
        outfile.write('countsTable <- subset(data, select=c('+', '.join([str(i) for i in range(5,5+len(BAM1)+len(BAM2))])+'))\n')
        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- c(' + ', '.join(['"'+LABEL1+'"']*len(BAM1)) + ', ' + ', '.join(['"'+LABEL2+'"']*len(BAM2)) + ')\n')
        outfile.write('cds <- newCountDataSet( countsTable, conds )\n')
        outfile.write('cds <- estimateSizeFactors( cds )\n')
        outfile.write('sizeFactors(cds)\n')                                                               
        outfile.write('cds <- estimateDispersions( cds ,method="blind", sharingMode="fit-only")\n') ##without replicates                        
        outfile.write('res <- nbinomTest( cds, "'+LABEL1+'", "'+LABEL2+'" )\n')
        outfile.write('rownames(res) <- res$id\n')                      
        outfile.write('write.table(res, file = "'+FILEDIR+'DESeq.res.txt", append = FALSE, sep= "\t" )\n')
        outfile.write('sink()')

def run():
    write_script()
    os.system("R < " + FILEDIR + "DESeq.R --no-save")
