**Project:** Gene expression and eQTL analysis reflect the heterogeneity in the inflammatory status of the duodenal epithelial lining in coeliac disease <br><br>
**Author:** Aarón D. Ramírez-Sánchez<br>
**PI:** Iris Jonkers<br>
Department of Genetics, University Medical Center Groningen, Netherlands<br>


This repository contains the scripts used for raw data processing, data preprocessing, and downstream analysis.

This file describes the pipeline followed for the Epithelial CeD project.

## Preparations:

00_Formatting_metadata.R: Creating a metadata file with the info necessary for downstream analysis and a file of reference to change sample names with new nomitation.

00_Formatting_count_matrix.R: Creating RNA count data with new names and filter out outliers

01_DifExpAnalysis_DESeq.R: Creating DSEQ object and differential expression analysis tables

02_Enrichment_analysis.R: Creating enrichment analysis tables and objects. Setting cutoffs for p vals and L2FC. The output of this script is used on 03

03_Enrichment_plots.R: Creating enrichment analysis dotplots

04_ZeroVar_count_matrix.R: Creating a matrix filtering low expression and zero variance genes

05_Heatmaps_and_clusters.R: Setting cutoffs for p vals and L2FC

06_Enrichment_analysis_clusters.R: Creating enrichment analysis tables and objects for clusters of genes. The outputs of this script are used in 07.

07_Enrichment_plots_clusters.R

08_PCA_plots.R

09_DifExpAnalysis_DESeq_clusters.R

10_Enrichment_analysis_clusters_DEG.R:

11_Enrichment_plots_clusters_DEG.R:

12_Heatmaps_and_clusters_clusters_DEG.R: Not completed, maybe not interesting

13_DEG_overview_plots.R

14_DEG_clusters_overview_plot.R
