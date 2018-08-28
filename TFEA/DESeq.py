__author__ = 'Jonathan Rubin'

import os, config, matplotlib.pyplot as plt

def write_script():
    if (len(config.BAM1) > 1 and len(config.BAM2) > 1):
        outfile = open(config.FILEDIR + 'DESeq.R','w')
        outfile.write('sink("'+config.FILEDIR+'DESeq.Rout")\n')
        outfile.write('library("DESeq2")\n')
        outfile.write('data <- read.delim("'+config.COUNT_FILE+'", sep="\t", header=TRUE)\n')
        outfile.write('countsTable <- subset(data, select=c('+', '.join([str(i) for i in range(5,5+len(config.BAM1)+len(config.BAM2))])+'))\n')
        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- as.data.frame(c(' + ', '.join(['"'+config.LABEL1+'"']*len(config.BAM1)) + ', ' + ', '.join(['"'+config.LABEL2+'"']*len(config.BAM2)) + '))\n')
        outfile.write('colnames(conds) <- c("treatment")\n')
        outfile.write('ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsTable, colData = conds, design = ~ treatment)\n')
        outfile.write('dds <- DESeq(ddsFullCountTable)\n')
        outfile.write('res1 <- results(dds,alpha = 0.05, contrast=c("treatment","'+config.LABEL2+'","'+config.LABEL1+'"))\n')
        outfile.write('resShrink <- lfcShrink(dds, res = res1, contrast = c("treatment","'+config.LABEL2+'","'+config.LABEL1+'"))\n')
        outfile.write('resShrink$fc <- 2^(resShrink$log2FoldChange)\n')
        outfile.write('res <- resShrink[c(1:3,7,4:6)]\n')
        outfile.write('write.table(res, file = "'+config.FILEDIR+'DESeq.res.txt", append = FALSE, sep= "\t" )\n')
        outfile.write('sink()')
    else:
        outfile = open(config.FILEDIR + 'DESeq.R','w')
        outfile.write('sink("'+config.FILEDIR+'DESeq.Rout")\n')
        outfile.write('library("DESeq")\n')
        outfile.write('data <- read.delim("'+config.COUNT_FILE+'", sep="\t", header=TRUE)\n')
        outfile.write('countsTable <- subset(data, select=c('+', '.join([str(i) for i in range(5,5+len(config.BAM1)+len(config.BAM2))])+'))\n')
        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- c(' + ', '.join(['"'+config.LABEL1+'"']*len(config.BAM1)) + ', ' + ', '.join(['"'+config.LABEL2+'"']*len(config.BAM2)) + ')\n')
        outfile.write('cds <- newCountDataSet( countsTable, conds )\n')
        outfile.write('cds <- estimateSizeFactors( cds )\n')
        outfile.write('sizeFactors(cds)\n')                                                               
        outfile.write('cds <- estimateDispersions( cds ,method="blind", sharingMode="fit-only")\n') ##without replicates                        
        outfile.write('res <- nbinomTest( cds, "'+config.LABEL1+'", "'+config.LABEL2+'" )\n')
        outfile.write('rownames(res) <- res$id\n')                      
        outfile.write('write.table(res, file = "'+config.FILEDIR+'DESeq.res.txt", append = FALSE, sep= "\t" )\n')
        outfile.write('sink()')

def run():
    write_script()
    os.system("R < " + config.FILEDIR + "DESeq.R --no-save")
    plot_MA(config.FILEDIR+'DESeq.res.txt')

def plot_MA(deseq_file):
    x = list()
    sigx = list()
    y = list()
    sigy = list()
    with open(deseq_file,'r') as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            basemean = float(line[2])
            try:
                log2fc = float(line[6])
                padj = float(line[-1])
                if padj < config.PADJCUTOFF:
                    sigx.append(basemean)
                    sigy.append(log2fc)
                else:
                    x.append(basemean)
                    y.append(log2fc)
            except:
                pass

    #Creates an MA-Plot of the region expression
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    ax.scatter(x=x,y=y,color='black',edgecolor='')
    ax.scatter(x=sigx,y=sigy,color='red',edgecolor='')
    ax.set_title("DE-Seq MA-Plot",fontsize=14)
    ax.set_ylabel("Log2 Fold-Change ("+config.LABEL1+"/"+config.LABEL2+")",fontsize=14)
    ax.set_xlabel("Average Expression",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    # plt.show()
    plt.savefig(config.FIGUREDIR + 'DESEQ_MA_Plot.png',bbox_inches='tight')

if __name__ == "__main__":
    deseq_file = '/Users/jonathanrubin/Google_Drive/Colorado University/Jonathan/TFEA_outputs/Allen2014/TFEA_DMSO1-DMSO2_3/temp_files/DESeq.res.txt'
    plot_MA(deseq_file)
    plt.cla()