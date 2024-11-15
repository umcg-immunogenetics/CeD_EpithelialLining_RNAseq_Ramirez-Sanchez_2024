#Script generated by: Aaron D. Ramirez-Sanchez 
#Goal: To format metadata and create a sample_sample conversion

##ENVIRONMENT
library(ggplot2)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(eulerr)
library(patchwork)

##FUNCTIONS
#####PLOTS

theme.plain <- function(p, base_size = 10, base_family = "ArialMT") {
  p <- p + theme_linedraw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          #panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", size=0.75),
          axis.ticks = element_line(size=0.75),
          axis.text = element_text(size=base_size, family="ArialMT", face="plain"),
          strip.background = element_blank(),
          strip.text = element_text(size=base_size, family="ArialMT", face="plain"),
          legend.key = element_blank(),
          legend.text = element_text(size=base_size, family = "ArialMT", face="plain"),
          complete = TRUE,
          plot.title = element_text(hjust=0.5))
  return(p)
}


##INPUT
coldata <- read.delim(file = "./data/processed/metadata/extended_metadata.tsv", sep = "\t", row.names = 1)

#reading merged DEG 
#Loading merged DEG 
DEG.merged <- read.csv(file = "./output/tables/DEG/comparisons_merged.csv", header = T, stringsAsFactors = F)
DEG.merged.clusters <- read.csv(file = "./output/tables/DEG/clusters_comparisons_merged.csv", header = T, stringsAsFactors = F)
DEG.merged <- rbind(DEG.merged, DEG.merged.clusters)

#Filters used for the DEG
L2FC.filter <- 1.0
p.adj.filter <- 0.01

RNA.ref <- read.table("/Users/R_projects/Oslo/results/eQTL/references/Homo_sapiens.GRCh37.75.genes.bed", sep = "\t", header = F)
names(RNA.ref) <- c("chr", "start_pos", "end_pos", "strand", "ENSEMBL", "SYMBOL")
RNA.ref <- RNA.ref[!duplicated(RNA.ref$ENSEMBL),]
Y.genes <- RNA.ref[RNA.ref$chr == "Y", ]$ENSEMBL

#colors: Danger Red, Tap Shoe and Blue Blossom
pallete_condition <- c("#D9514E", "#2A2B2D", "#2DA8D8")
#colors: Danger Red and 1/2, Tap Shoe and 1/2; and Blue Blossom and 1/2
pallete_condition_updown <- c("#D9514E","#FFACAA", "#2A2B2D",
                              "#6B6E73", "#2DA8D8","#7ECDEA", 
                              "#2A2B2DFF", "#8B8090", "#315261",
                              "#2DA8D8FF", "#00D5C9", "#5DE7AA",
                              "#D9514EFF", "#BC519B", "#9260B5", 
                              "#DC8AAA")

##MAIN

#Filtering data
DEG.merged <- DEG.merged[(abs(DEG.merged$log2FoldChange) > L2FC.filter & !is.na(DEG.merged$log2FoldChange)) & 
                           (DEG.merged$padj < p.adj.filter & !is.na(DEG.merged$padj)), ]
DEG.merged$Direction <- factor(DEG.merged$Direction, levels = c("Upregulated", "Downregulated"))
DEG.merged <- DEG.merged[!DEG.merged$ENSEMBL %in% Y.genes, ]

#plot of genes up and down
p <- ggplot(DEG.merged, aes(x = Comparison, fill = Direction)) +
  geom_bar(position = "dodge") +
  geom_text(aes(label = ..count..), stat = "count", position = position_dodge(0.9), vjust = -0.5) +
  scale_fill_manual(values = c("#b2182b","#2166ac"))

#Euler diagram
#Creating function to create matrix with the DEG 
TCDvsCTRL <- DEG.merged[DEG.merged$Comparison == "TCDvsCTRL",]
UCDvsCTRL <- DEG.merged[DEG.merged$Comparison == "UCDvsCTRL",]
UCDvsTCD <- DEG.merged[DEG.merged$Comparison == "UCDvsTCD",]
G2.UCDvsTCD <- DEG.merged[DEG.merged$Comparison == "G2.UCDvsTCD",]
UCD.3vs2 <- DEG.merged[DEG.merged$Comparison == "UCD.3vs2",]
TCD.2vs1 <- DEG.merged[DEG.merged$Comparison == "TCD.2vs1",]

list.up <- list(TCDvsCTRL = TCDvsCTRL[TCDvsCTRL$Direction == "Upregulated",]$ENSEMBL,
                UCDvsCTRL = UCDvsCTRL[UCDvsCTRL$Direction == "Upregulated",]$ENSEMBL,
                UCDvsTCD = UCDvsTCD[UCDvsTCD$Direction == "Upregulated",]$ENSEMBL,
                G2.UCDvsTCD = G2.UCDvsTCD[G2.UCDvsTCD$Direction == "Upregulated",]$ENSEMBL,
                UCD.3vs2 = UCD.3vs2[UCD.3vs2$Direction == "Upregulated",]$ENSEMBL,
                TCD.2vs1 = TCD.2vs1[TCD.2vs1$Direction == "Upregulated",]$ENSEMBL
)

list.down <- list(TCDvsCTRL = TCDvsCTRL[TCDvsCTRL$Direction == "Downregulated",]$ENSEMBL,
                  UCDvsCTRL = UCDvsCTRL[UCDvsCTRL$Direction == "Downregulated",]$ENSEMBL,
                  UCDvsTCD = UCDvsTCD[UCDvsTCD$Direction == "Downregulated",]$ENSEMBL,
                  G2.UCDvsTCD = G2.UCDvsTCD[G2.UCDvsTCD$Direction == "Downregulated",]$ENSEMBL,
                  UCD.3vs2 = UCD.3vs2[UCD.3vs2$Direction == "Downregulated",]$ENSEMBL,
                  TCD.2vs1 = TCD.2vs1[TCD.2vs1$Direction == "Downregulated",]$ENSEMBL
)

m.up.genes <- list_to_matrix(list.up)

m.down.genes <- list_to_matrix(list.down)


fit.up <- euler(m.up.genes, shape = "ellipse")
fit.down <- euler(m.down.genes, shape = "ellipse")

colnames(m.up.genes)
rownames(fit.up)

p.eu.up <- plot(fit.up,
                quantities = TRUE,
                edges = F,
                fills = pallete_condition_updown[c(5,1,2, 14, 15, 11)],
                labels = list(font = 2, fontsize = 12),
                main = "Upregulated")
p.eu.up.bw <- plot(fit.up,
                quantities = TRUE,
                edges = T,
                fills = F,
                labels = list(font = 2, fontsize = 12),
                main = "Upregulated")

p.eu.down <- plot(fit.down,
                  quantities = TRUE,
                  edges = F,
                  fills = pallete_condition_updown[c(5,1,2, 14, 15, 11)],
                  labels = list(font = 2, fontsize = 12),
                  main = "Downregulated")

p.eu.down.bw <- plot(fit.down,
                  quantities = TRUE,
                  edges = T,
                  fills = F,
                  labels = list(font = 2, fontsize = 12),
                  main = "Downregulated")
c("#D9514E", "#FFACAA", "#2A2B2D", "#6B6E73", "#2DA8D8", "#7ECDEA")

##UPSET plots
# compare two UpSet plots

m.comb.up <-  make_comb_mat(list.up)
m.comb.down <-  make_comb_mat(list.down)
set_name(m.comb.up)
comb_name(m.comb.up)
Unique.up.UCD3vs2 <- extract_comb(m.comb.up, "000010")
Unique.down.UCD3vs2 <- extract_comb(m.comb.down, "000010")

DEG.unique.UCD3vs2 <- rbind(DEG.merged[DEG.merged$ENSEMBL %in% Unique.up.UCD3vs2, ], DEG.merged[DEG.merged$ENSEMBL %in% Unique.down.UCD3vs2, ])

max1 <- max(c(set_size(m.comb.up), set_size(m.comb.down)))
max2 <- max(c(comb_size(m.comb.up), comb_size(m.comb.down)))

upset.plot <- UpSet(m.comb.up, 
                    top_annotation = upset_top_annotation(m.comb.up, ylim = c(0, max2)),
                    right_annotation = upset_right_annotation(m.comb.up, ylim = c(0, max1)),
                    column_title = "Upregulated") +
  UpSet(m.comb.down, 
        top_annotation = upset_top_annotation(m.comb.down, ylim = c(0, max2)),
        right_annotation = upset_right_annotation(m.comb.down, ylim = c(0, max1)),
        column_title = "Downregulated")



##OUTPUT

pdf(file = "./output/plots/DEG_clusters_and_general_overview.pdf", width=6, height=4.5,  family = "ArialMT")
theme.plain(p)+theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf(file = "./output/plots/DEG_clusters_and_general_up_euler.pdf", width=6, height=6,  family = "ArialMT")
p.eu.up
p.eu.up.bw
dev.off()

pdf(file = "./output/plots/DEG_clusters_and_general_down_euler.pdf", width=6, height=6,  family = "ArialMT")
p.eu.down
p.eu.down.bw
dev.off()

pdf(file = "./output/plots/DEG_clusters_and_general_upset.pdf", width=9, height=4,  family = "ArialMT")
upset.plot
dev.off()

write.table(DEG.unique.UCD3vs2 , file = "./output/tables/DEG/DEG.unique.UCD3vs2.txt",
            row.names = F, quote = F, sep = "\t")





# ERICHMENT ANALYSIS

#REACTOME enrichment
DEG.unique.UCD3vs2$Comparison <- gsub("\\.", "_", DEG.unique.UCD3vs2$Comparison)

DEG.Pathway <- compareCluster(ENTREZID~Direction, data=DEG.unique.UCD3vs2, fun="enrichPathway", readable=T)
DEG.Pathway@compareClusterResult$minuslog10p.adjust <- -log10(DEG.Pathway@compareClusterResult$p.adjust)


#GO: Molecular function
DEG.mf <- compareCluster(ENTREZID~Direction, data=DEG.unique.UCD3vs2, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "MF", readable=T) 
DEG.mf@compareClusterResult$minuslog10p.adjust <- -log10(DEG.mf@compareClusterResult$p.adjust)


#GO: Biological process
DEG.bp <- compareCluster(ENTREZID~Direction, data=DEG.unique.UCD3vs2, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", readable=T) 
DEG.bp@compareClusterResult$minuslog10p.adjust <- -log10(DEG.bp@compareClusterResult$p.adjust)


#GO: Celullar component
DEG.cc <- compareCluster(ENTREZID~Direction, data=DEG.unique.UCD3vs2, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "CC", readable=T) 
DEG.cc@compareClusterResult$minuslog10p.adjust <- -log10(DEG.cc@compareClusterResult$p.adjust)


#KEGG enrichment
DEG.kegg <- compareCluster(ENTREZID~Direction, data=DEG.unique.UCD3vs2, fun="enrichKEGG") 
DEG.kegg@compareClusterResult$minuslog10p.adjust <- -log10(DEG.kegg@compareClusterResult$p.adjust)

p1 <- dotplot(DEG.Pathway, 
        showCategory=5, 
        color = "minuslog10p.adjust", 
        font.size = 9, label_format = 50) + 
  scale_colour_gradient(low="#D0E0F3", high="#257EE4") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- dotplot(DEG.mf, 
        showCategory=5, 
        color = "minuslog10p.adjust", 
        font.size = 9, label_format = 50) + 
  scale_colour_gradient(low="#D0E0F3", high="#257EE4") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- dotplot(DEG.bp, 
        showCategory=5, 
        color = "minuslog10p.adjust", 
        font.size = 9, label_format = 50) + 
  scale_colour_gradient(low="#D0E0F3", high="#257EE4") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 <- dotplot(DEG.cc, 
        showCategory=5, 
        color = "minuslog10p.adjust", 
        font.size = 9, label_format = 50) + 
  scale_colour_gradient(low="#D0E0F3", high="#257EE4") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(file = "./output/plots/DEG_unique_UCDg3vs2_enrichments.pdf", width=9, height=4,  family = "ArialMT")
theme.plain(p1) +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(size = 8, lineheight = 0.7), 
        strip.text.x = element_text(size = 9, vjust=0.5)) +
  labs(colour = "-log10(adj. P-value)") +
  ggtitle("Pathway Enrichment")

theme.plain(p2) +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(size = 8, lineheight = 0.7), 
        strip.text.x = element_text(size = 9, vjust=0.5)) +
  labs(colour = "-log10(adj. P-value)") +
  ggtitle("GO: Molecular function Enrichment")

theme.plain(p3) +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(size = 8, lineheight = 0.7), 
        strip.text.x = element_text(size = 9, vjust=0.5)) +
  labs(colour = "-log10(adj. P-value)") +
  ggtitle("GO: Biological process Enrichment")

theme.plain(p4) +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(size = 8, lineheight = 0.7), 
        strip.text.x = element_text(size = 9, vjust=0.5)) +
  labs(colour = "-log10(adj. P-value)") +
  ggtitle("GO: Celullar component Enrichment")
dev.off()