decon_eQTLs <- read.table("./output/tables/decon-eQTLs/decon_FACS2_RNA_marsh.txt", header = T)
decon_eQTLs <- decon_eQTLs[decon_eQTLs$sig.Epithelial | decon_eQTLs$sig.Immune,]

adding.sig.eQTL.for2 <- function(df) {
  df$sig.Epithelial <- df$Epithelial_cells_pvalue < 0.01
  df$sig.Immune <- df$Immune_Cells_pvalue < 0.01
  return(df)
}

decon_eQTLs_01 <- adding.sig.eQTL.for2(decon_eQTLs)
decon_eQTLs_01 <- decon_eQTLs_01[decon_eQTLs_01$sig.Epithelial | decon_eQTLs_01$sig.Immune,]

SNPs_01 <- decon_eQTLs_01[!duplicated(decon_eQTLs_01$SNP),]$SNP
Genes_01 <- decon_eQTLs_01[!duplicated(decon_eQTLs_01$SYMBOL),]$SYMBOL

write.table(SNPs_01, 
            file="./output/tables/decon-eQTLs/SNPs_mss_01.txt", 
            row.names = F, quote = F, col.names =  F, sep = "\t")

write.table(Genes_01, 
            file="./output/tables/decon-eQTLs/Genes_mss_01.txt", 
            row.names = F, quote = F, col.names =  F, sep = "\t")

write.table(decon_eQTLs_01, 
            file="./output/tables/decon-eQTLs/eQTLs_mss_01.txt", 
            row.names = F, quote = F, sep = "\t")

#for ssh 
head -n1 SNPs_mss_01.txt genotypes_decon-eQTL.txt > genotypes_decon-eQTL_toplot1.txt
awk -F "\t" 'FNR==NR { a[$1]; next } ($1 in a)' SNPs_mss_01.txt genotypes_decon-eQTL.txt >> genotypes_decon-eQTL_toplot1.txt
awk -F "\t" '{$1=""; print $0}' genotypes_decon-eQTL_toplot1.txt > genotypes_decon-eQTL_toplot2.txt