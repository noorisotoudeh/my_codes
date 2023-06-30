# Qiagen
This repository provides example codes used for some basic and advanced data analysis in ...
# 1- single cell rna-seq
in this directory you can find the notebooks for data analysis used for single cell rna-seq data in [ETP-ALL](https://ashpublications.org/blood/article/137/18/2463/474247/Single-cell-RNA-seq-reveals-developmental) or [myeloma](https://www.nature.com/articles/s41556-021-00766-y) paper. Note that, since these analysis has been done when popular package Seurat was not as developed as what it is now, the functions like SCTransform didn't use and instead SingleCellExperiment objects and its related packages applied for the most of analysis.
# 2- bulk rna-seq
Deseq2 and edgeR are two well known packages that provides methods to test for differential gene expressions. The notebook contains DGE analysis using Deseq2 package for bulk rna samples.
# 3- milti omics
here I will show how to use combined scRNA-seq and single cell chromatin accessibility (scATAC-seq) data to create gene regulatory networks using scenic+ package. The final outputs can be displayed in cytoscape (desktop/web version).
# 4- shiny app
this is a development of ShinyCell package and it gives you an example in spatial ATAC-seq data analysis. 
