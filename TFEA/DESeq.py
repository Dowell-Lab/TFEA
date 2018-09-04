__author__ = 'Jonathan Rubin'

import os, math, config, matplotlib.pyplot as plt, numpy as np

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

    x = [x for _,x in sorted(zip(up_p,up_x))] + [x for _,x in sorted(zip(dn_p,dn_x), reverse=True)]
    y = [y for _,y in sorted(zip(up_p,up_y))] + [y for _,y in sorted(zip(dn_p,dn_y), reverse=True)]
    c = plt.cm.RdYlGn(np.linspace(0,1,len(x)))
    # c =np.linspace(0,1,len(x))
    # cm = plt.cm.get_cmap('RdYlGn')
    #Creates an MA-Plot of the region expression
    F = plt.figure(figsize=(7,6))
    ax = plt.subplot(111)
    plt.scatter(x=x,y=y,color=c,edgecolor='')
    ax.set_title("DE-Seq MA-Plot",fontsize=14)
    ax.set_ylabel("Log2 Fold-Change ("+config.LABEL2+"/"+config.LABEL1+")",fontsize=14)
    ax.set_xlabel("Log10 Average Expression",fontsize=14)
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='on')
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
    # plt.colorbar()
    # plt.show()
    plt.savefig(config.FIGUREDIR + 'DESEQ_MA_Plot.png',bbox_inches='tight')

if __name__ == "__main__":
    deseq_file = '/Users/jonathanrubin/Google_Drive/Colorado University/Jonathan/TFEA_outputs/Allen2014/TFEA_DMSO-NUTLIN_3/temp_files/DESeq.res.txt'
    plot_MA(deseq_file)
    plt.cla()