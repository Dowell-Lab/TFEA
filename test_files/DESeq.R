library("DESeq2")
data <- read.delim("/Users/joru1876/Google_Drive/Colorado_University/Jonathan/TFEA/test_files/count_file.header.bed", sep="	", header=TRUE)
countsTable <- subset(data, select=c(5, 6, 7, 8))

rownames(countsTable) <- data$region
conds <- as.data.frame(c("condition1", "condition1", "condition2", "condition2"))

colnames(conds) <- c("treatment")
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countsTable,
                                            colData = conds,
                                            design = ~ treatment)

dds <- DESeq(ddsFullCountTable)
res <- results(dds, alpha = 0.05, contrast=c("treatment", "condition2",
                                                            "condition1"))
res$fc <- 2^(res$log2FoldChange)
res <- res[c(1:3,7,4:6)]

write.table(res, file = "/Users/joru1876/Google_Drive/Colorado_University/Jonathan/TFEA/test_files/DESeq.res.txt", append = FALSE, sep= "	" )
sink()