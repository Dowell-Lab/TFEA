library("DESeq")
data <- read.delim("/Users/joru1876/Google_Drive/Colorado_University/Jonathan/TFEA/TFEA/test/test_files/test_output/temp_files/count_file.header.bed", sep="	", header=TRUE)

countsTable <- subset(data, select=c(5, 6))

rownames(countsTable) <- data$region
conds <- c("DMSO", "Nutlin")

cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
sizeFactors(cds)                                                               
cds <- estimateDispersions( cds ,method="blind",                         sharingMode="fit-only")

res <- nbinomTest( cds, "DMSO", "Nutlin" )
rownames(res) <- res$id                      
write.table(res, file = "/Users/joru1876/Google_Drive/Colorado_University/Jonathan/TFEA/TFEA/test/test_files/test_output/temp_files/DESeq.res.txt", append = FALSE, sep= "	" )