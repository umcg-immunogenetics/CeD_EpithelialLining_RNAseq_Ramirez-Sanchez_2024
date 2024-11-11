CeD.snps <- read.delim(file = "./data/raw/reference/Top_CeD_SNPs_May2019.txt", sep = "\t", header = T)
CeD.snps$position.2 <- gsub("chr","",CeD.snps$hg19_position)
CeD.snps$chr <- gsub(":.*","",CeD.snps$position.2)
CeD.snps$pos1 <- as.integer(gsub(".*:","",CeD.snps$position.2))
CeD.snps$pos2 <- CeD.snps$pos1+1
names(CeD.snps)
CeD.bed <- CeD.snps[,c("chr", "pos1", "pos2", "Top_SNP", "MAF_1000G_Eur", "Locus")]
names(CeD.bed) <- NULL
write.table(CeD.bed, 
            './data/processed/reference/Top_CeD_SNPs_May2019.bed', 
            sep='\t', row.names = F, quote=F)
