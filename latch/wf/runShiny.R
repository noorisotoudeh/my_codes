library("ArchR")
library("GenomicRanges")
library('BSgenome')
library("dplyr")                                    # Load dplyr package
library("plyr")                                     # Load plyr package
library("readr")
library(qdap)
library('Seurat')
library(ShinyCell)
library('BSgenome.Mmusculus.UCSC.mm10')
library('BSgenome.Hsapiens.UCSC.hg38')
suppressPackageStartupMessages(require(tidyverse));



dir.create("/root/results", showWarnings = FALSE)
setwd("/root/results")

rawPath <- "/root/"
dataFiles <- dir(rawPath, "*.R$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)




find_func <- function(tempdir,pattern){
    
   list.files(
  path = tempdir, # replace with the directory you want
  pattern = pattern, # has "test", followed by 0 or more characters,
                             # then ".csv", and then nothing else ($)
  full.names = TRUE # include the directory in the result
        , recursive = TRUE
)
    
}


args <- commandArgs(trailingOnly=TRUE)
tempdir <- args[1]

ArchRobj <- system(paste0("find ", tempdir, " -name '*_ArchRProject' -type d"), intern = TRUE)


proj3 <- loadArchRProject(path = ArchRobj, force = FALSE, showLogo = TRUE)

############------------------------Identifying Marker Genes--------------------------###################

markersGS <- getMarkerFeatures(
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",#"wilcoxon"
  
)

saveRDS(markersGS,"markersGS_clusters.rds")

for (rd in names(proj3@reducedDims)){
  if (rd=="Harmony"){
    proj3 <- addImputeWeights(proj3,reducedDims = "Harmony")
    
  } else {
    proj3 <- addImputeWeights(proj3)
    
  }
}
# save for shiny

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20", plotLog2FC = TRUE,
  #   labelMarkers = seleceted_markers,
  transpose = F,  returnMatrix = TRUE
)

write.csv(heatmapGS,"genes_per_cluster_hm.csv")


# per sample

######------------------------------Identifying Marker Genes--------------------------------#######

markersGS <- getMarkerFeatures(
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",#"wilcoxon"
)

# save for shiny app
saveRDS(markersGS,"markersGS_sample.rds")


for (rd in names(proj3@reducedDims)){
  if (rd=="Harmony"){
    proj3 <- addImputeWeights(proj3,reducedDims = "Harmony")
    
  } else {
    proj3 <- addImputeWeights(proj3)
    
  }
}


# save for shiny

if (length(unique(proj3$Sample))>1){
  
  # save for shiny
  
  heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS, 
    cutOff = "FDR <= 0.05 & Log2FC >= 0.20", plotLog2FC = TRUE,
    #   labelMarkers = seleceted_markers,
    transpose = F,  returnMatrix = TRUE
  )
  
} else {
  
  heatmapGS <- "there is not enough samples to be compared with!"  
}   

write.csv(heatmapGS,"genes_per_sample_hm.csv")

# per treatment

######------------------------------Identifying Marker Genes--------------------------------#######

markersGS <- getMarkerFeatures(
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",#"wilcoxon"
)
# save for shiny app
saveRDS(markersGS,"markersGS_treatment.rds")


for (rd in names(proj3@reducedDims)){
  if (rd=="Harmony"){
    proj3 <- addImputeWeights(proj3,reducedDims = "Harmony")
    
  } else {
    proj3 <- addImputeWeights(proj3)
    
  }
}


# save for shiny

if (length(unique(proj3$Condition))>1){
  
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20", plotLog2FC = TRUE,
  #   labelMarkers = seleceted_markers,
  transpose = F,  returnMatrix = TRUE
)

} else {
  
  heatmapGS <- "there is not enough Conditions to be compared with!"  
}   


write.csv(heatmapGS,"genes_per_treatment_hm.csv")                     



# Volcano plots for genes
if (length(unique(proj3$Condition))>1){
  
ncells <- length(proj3$cellNames)

markerList <- getMarkerFeatures(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),maxCells = ncells ,normBy = "none",
  testMethod = "ttest")

req_DF <- as.data.frame(getCellColData(proj3))
df1 <- table(req_DF$Clusters,req_DF$Condition)
distr <- as.data.frame.matrix(round(prop.table(as.matrix(df1),1),2))
lst <- list()

for(i in 1:nrow(distr)) {
  row <- distr[i,]
  if (
    sum(unname(unlist(row))>= 0.90) == 1) {
    rownames(row) -> lst[[i]]
  }
}
not_req_list <- unlist(lst)

req_clusters <- unique(proj3$Clusters)
req_clusters <- req_clusters[order(as.numeric(gsub("C","",req_clusters)))]
req_clusters <- req_clusters[which(!req_clusters%in%not_req_list)]

markerList_C <- list()
proj_C <- list()

for (i in seq_along(req_clusters)) {
  
  idxSample <- BiocGenerics::which(proj3$Clusters == req_clusters[i])
  
  cellsSample <- proj3$cellNames[idxSample]
  proj_C[i] <- proj3[cellsSample,]
  
  ncells[i] <- length(proj_C[[i]]$cellNames)
  
  # per each cluster separately
  markerList_C[[i]] <- getMarkerFeatures(
    ArchRProj = proj_C[[i]],
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Condition",
    bias = c("TSSEnrichment", "log10(nFrags)"),maxCells = ncells[[i]] ,normBy = "none",
    testMethod = "ttest")
}
names(markerList_C) <- req_clusters


gsm <- getMatrixFromProject(proj3)
gsm_mat <- assay(getMatrixFromProject(proj3),"GeneScoreMatrix")

which(rowSums(is.na(gsm_mat))>0)
any(rowSums((gsm_mat))==0)

empty_gene_idx <- which(rowSums((gsm_mat))==0)
empty_gene <- rowData(gsm)$name[empty_gene_idx]



markerList_df1 <- assay(markerList, "Log2FC")
markerList_df2 <- assay(markerList, "Pval")
markerList_df3 <- assay(markerList, "FDR")
markerList_df <- cbind(markerList_df1,markerList_df2,markerList_df3)
markerList_df$genes<- rowData(markerList)$name
markerList_df$cluster <- rep("All",length(rownames(markerList_df)))

# # we only want to see results of one set , say just sham
markerList_df <- markerList_df[,c(1,3,5,7,8)]
colnames(markerList_df) <- c("avg_log2FC","p_val","p_val_adj","gene","cluster")

# percluster


markerList_df1_C <- list()
markerList_df2_C <- list()
markerList_df3_C <- list()
markerList_df_C <- list()

# for (i in (1:nClust)){
for (i in seq_along(req_clusters)){
  
  cluster <- req_clusters[i]
  
  markerList_df1_C[[i]] <- assay(markerList_C[[i]], "Log2FC")
  markerList_df2_C[[i]] <- assay(markerList_C[[i]], "Pval")
  markerList_df3_C[[i]] <- assay(markerList_C[[i]], "FDR")
  
  markerList_df_C[[i]] <- cbind(markerList_df1_C[[i]]
                                ,markerList_df2_C[[i]]
                                ,markerList_df3_C[[i]])
  
  markerList_df_C[[i]]$genes<- rowData(markerList_C[[i]])$name
  markerList_df_C[[i]]$cluster <- rep(cluster,length(rownames(markerList_df_C[[i]])))
  
  # # we only want to see results of one set , say just sham
  markerList_df_C[[i]] <- markerList_df_C[[i]][,c(1,3,5,7,8)]
  colnames(markerList_df_C[[i]]) <- c("avg_log2FC","p_val","p_val_adj","gene","cluster")
}


names(markerList_df_C) <- req_clusters

markersGS_merged_df <- do.call("rbind", markerList_df_C)

# also data frame for all clusters together needs to be added

markersGS_merged_df <- rbind(markerList_df,markersGS_merged_df)

# remove empty genes
markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$gene%in%empty_gene),]

# remove na values
markersGS_merged_df <- na.omit(markersGS_merged_df)

# remove FDR equal to 0
markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$p_val_adj== 0),]



# make logfc limiation between 1 and -1

markersGS_merged_df <- markersGS_merged_df[which(abs(markersGS_merged_df$avg_log2FC)< 1.2),]

markersGS_merged_df$Significance = ifelse(markersGS_merged_df$p_val < 10^-2 , 
                                          #                                           abs(markersGS_merged_df$avg_log2FC) >= 0.58, 
                                          ifelse(markersGS_merged_df$avg_log2FC> 0.0 
                                                 ,colnames(markerList)[1],colnames(markerList)[2]),
                                          'Not siginficant')

de <- markersGS_merged_df

} else {
  
  de <- "there is not enough samples to be compared with!"  
}   


write.table(de,"inpMarkers.txt", sep = '\t', quote = F, row.names = F)


suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library(data.table))

tempdir <- "/root/results"
genes_per_cluster_hm <- find_func(tempdir,"genes_per_cluster_hm.csv")
hm_per_clust <- read.csv(genes_per_cluster_hm)

genes_per_sample_hm <- find_func(tempdir,"genes_per_sample_hm.csv")
hm_per_sample <- read.csv(genes_per_sample_hm)

genes_per_cond_hm <- find_func(tempdir,"genes_per_treatment_hm.csv")
hm_per_cond <- read.csv(genes_per_cond_hm)



nClust = length(unique(proj3$Clusters))


df = list()


for (i in seq_along(1:nClust)){
  df[[i]] <- hm_per_clust[,c(1,i+1)]
  
  #select top 5 values by group
  
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:5,1]
  
}
final <- do.call(rbind, df)

req_genes1 <- unlist(df)

req_genes1<- req_genes1[!duplicated(req_genes1)]

# save the genes for default values in shiny app

write.csv(req_genes1,"req_genes1.csv")



# per treatment
if (length(unique(proj3$Condition))>1){
  
nConds = 2


df = list()

for (i in seq_along(1:nConds)){
  df[[i]] <- hm_per_cond[,c(1,i+1)]
  
  #select top 20 values by group
  
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:20,1]
  
}
final <- do.call(rbind, df)
# final
req_genes2 <- unlist(df)

req_genes2<- req_genes2[!duplicated(req_genes2)]

} else {
  
  req_genes2 <- "there is not enough condition to be compared with!"  
}   


write.csv(req_genes2,"req_genes2.csv")


# per sample

if (length(unique(proj3$Sample))>1){
  

nSamples = length(unique(proj3$Sample))


df = list()

for (i in seq_along(1:nSamples)){
  df[[i]] <- hm_per_sample[,c(1,i+1)]
  
  #select top 10 values by group
  
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:10,1]
  
}
final <- do.call(rbind, df)
# final
req_genes3 <- unlist(df)

req_genes3<- req_genes3[!duplicated(req_genes3)]
} else {
  
  req_genes3 <- "there is not enough samples to be compared with!"  
}


write.csv(req_genes3,"req_genes3.csv")


UMAPHarmony <-getEmbedding(ArchRProj = proj3, embedding = "UMAP", returnDF = TRUE)
write.csv(UMAPHarmony,"UMAPHarmony.csv")



############------------------------Peak Calling using MACS2--------------------------###################
pathToMacs2 <- findMacs2()
proj3 <- addGroupCoverages(ArchRProj = proj3
                           , groupBy = "Clusters"
                           ,maxCells = 1500
                           ,force = T)

species <- getGenome(ArchRProj = proj3)
if (species=="BSgenome.Hsapiens.UCSC.hg38"){
  data_species <- "hg38"
} else if (species=="BSgenome.Mmusculus.UCSC.mm10") {
  
  data_species <- "mm10"
}



if (data_species == "hg38") {
  addArchRGenome(data_species)
  geneAnnotation <- getGeneAnnotation()
  genomeAnnotation <- getGenomeAnnotation()
  genomeSize = 3.3e+09
} else if (data_species == "mm10") {
  addArchRGenome(data_species)
  geneAnnotation <- getGeneAnnotation()
  genomeAnnotation <- getGenomeAnnotation()
  genomeSize = 3.0e+09
}

proj3 <- addReproduciblePeakSet(
  ArchRProj = proj3,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  genomeSize = genomeSize,maxPeaks = 300000,
  force = TRUE
  
)
proj3 <- addPeakMatrix(proj3,force = TRUE)
# Motif enrichment (Deviation)
if("Motif" %ni% names(proj3@peakAnnotation)){
  if (data_species == "BSgenome.Hsapiens.UCSC.hg38" || data_species == "BSgenome.Mmusculus.UCSC.mm10") {
    proj3 <- addMotifAnnotations(ArchRProj = proj3, motifSet = "cisbp", name = "Motif", force = TRUE)
  } else {
    proj3 <- addMotifAnnotations(ArchRProj = proj3, motifSet = "encode", name = "Motif", force = TRUE
                                 , species = getGenome(ArchRProj = proj3))
  }
}

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj3, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj3,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
)

motif_lst <- unique(rownames(enrichMotifs))
split_string <- strsplit(motif_lst, split = "\\(")
fun1 <- function(list, nth){
  sapply(list, `[` , 1)
}
req_motifs1 <- gsub("_","-",fun1(split_string))
req_motifs1 <- gsub(" ","",req_motifs1)

rownames(enrichMotifs) <- req_motifs1

saveRDS(enrichMotifs,"enrichMotifs_clusters.rds")

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 50, transpose = F,returnMatrix = TRUE)

motif_lst <- unique(rownames(heatmapEM))
split_string <- strsplit(motif_lst, split = "\\(")
fun1 <- function(list, nth){
  sapply(list, `[` , 1)
}
req_motifs1 <- gsub("_","-",fun1(split_string))
req_motifs1 <- gsub(" ","",req_motifs1)

rownames(heatmapEM) <- req_motifs1

write.csv(heatmapEM,"motif_per_cluster_hm.csv")

nClust = length(unique(proj3$Clusters))

df = list()

for (i in seq_along(1:nClust)){
  df[[i]] <- hm_per_clust[,c(1,i+1)]
  
  #select top 5 values by group
  
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:5,1]
  
}
final <- do.call(rbind, df)

req_motifs1 <- unlist(df)

req_motifs1<- req_motifs1[!duplicated(req_motifs1)]

# save the motifs for default values in shiny app

write.csv(req_motifs1,"req_motifs1.csv")

# Add Motifs Matrix and Projections

proj3 <- addBgdPeaks(proj3, force = TRUE)
proj3 <- addDeviationsMatrix(
  ArchRProj = proj3,
  peakAnnotation = "Motif",
  force = TRUE
)
markersMotifs <- getMarkerFeatures(
  ArchRProj = proj3,
  useMatrix = "MotifMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z"
)
# creat deviation score 

markerMotifsList <- getMarkers(markersMotifs, cutOff = "FDR < 0.9 & MeanDiff >= 0")

motifs <- list()
for (i in seq_len(length(markerMotifsList))) {
  if (length(markerMotifsList[[i]]$name)>1) {
    motifs <- c(motifs, markerMotifsList[[i]]$name)
    #     motifs <- c(motifs, c(cluster1112_mm$name,cluster10_mm$name))
  }
}

if (length(motifs)>1) {
  motifs <- unlist(motifs)
  motifs <- paste0('z:', motifs)
  motifs <- unique(motifs)
  
  
  proj3 <- addImputeWeights(proj3)
  
  dev_score <- getDeviation_ArchR(ArchRProj = proj3, name = motifs
                                  , imputeWeights = getImputeWeights(proj3))
  
  #   dev_score <- as.data.frame(assay(MotifMatrix,"z"))
  
  dev_score[is.na(dev_score)] <- 0 #min(dev_score, na.rm = TRUE)
  
}
dev_score2 <- dev_score[!is.infinite(rowSums(dev_score)),]
colnames(dev_score2) <- rownames(getCellColData(proj3))
# remove 0 deviations per All samples 
all_zero <- names(which(rowSums(dev_score2)==0))
dev_score2 <- dev_score2[which(!rownames(dev_score2)%in%c(all_zero)),]

# for creating motif object go to after reading 




# Motif Logo

library("seqLogo")
require(ggseqlogo)
library(ArchR)
library(chromVARmotifs)

data("human_pwms_v1")

PWMs <- getPeakAnnotation(proj3, "Motif")$motifs

PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range


PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range


saveRDS(ProbMatrices,"seqlogo.rds")

#############################################add genome tracks

# gsm <- getMatrixFromProject(proj3)
# markerGenes  <- rowData(gsm)$name

tempdir <- "/root/results"
markerGenes <- unique(c(
                    read.csv(find_func(tempdir,"req_genes1.csv"))$x,
                    read.csv(find_func(tempdir,"req_genes2.csv"))$x,
                    read.csv(find_func(tempdir,"req_genes3.csv"))$x
))


# markerGenes <- unique(c(
#                      read.csv(find_func(tempdir,"genes_per_cluster_hm.csv"))$X,
#                      read.csv(find_func(tempdir,"genes_per_sample_hm.csv"))$X,
#                      read.csv(find_func(tempdir,"genes_per_treatment_hm.csv"))$X
# ))

               

 p <- plotBrowserTrack(
    ArchRProj = proj3,
    groupBy = "Clusters",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000
 )

 saveRDS(p,"ptracks.rds")
 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save ArchR object
 
saveArchRProject(ArchRProj = proj3, outputDirectory = "/root/results/ArchRProject")
 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tempdir <- args[1]

runs <- find_func(tempdir, "_SeuratObj_.*\\.rds$")

runs_spl <- strsplit(runs,"\\/|_")
runs_spl

n <- length(runs_spl[[1]]) - 3

fun1 <-  function(lst,nth){
    sapply(lst, `[`, nth)

}

names <- fun1(runs_spl,n+1)

names(runs) <- names

all <-  list()
for (i in seq_along(runs)){
        all[[i]] <- readRDS(runs[[i]])
}

for (i in seq_along(names)){
    
 all[[i]] <- RenameCells(all[[i]], new.names = paste0(unique(all[[i]]@meta.data$Sample)
                                               ,"#",colnames(all[[i]]),"-1"
                                                        ))
 
}


inputs1 <- all
#inputs2 <- args[1]

tempdir <- "/root/results"
inputs3 <- find_func(tempdir,"UMAPHarmony.csv")
req_genes1 <- find_func(tempdir,"req_genes1.csv")
req_genes2 <- find_func(tempdir, "req_genes2.csv")
req_genes3 <- find_func(tempdir,"req_genes3.csv")
#req_motifs1 <- args[6]
#req_motifs2 <- args[7]
#req_motifs3 <- args[8]

project <- args[2]
groupBy <- args[3]



#===========================#===========================#===========================#=================

main_func <- function(seurat_lst){
find_samples_name <- function(seurat_lst){
    
    sapply(seq_along(seurat_lst), function(i) unique(seurat_lst[[i]]@meta.data$Sample))

    
}

samples <- find_samples_name(seurat_lst)
    
D00_fun <- function(seurat_lst){
   toRemove <- lapply(seurat_lst, function(x) {names(which(colSums(is.na(x@assays[[1]]@counts))>0))}) 
    mapply(function(x,y) x[,!colnames(x) %in% y],seurat_lst,toRemove)
}


D00 <- D00_fun(seurat_lst)
Spatial_D00_fun <- function(D00){
    
Spatial_D00 <- lapply(D00, function(x) as.data.frame(x@images[[1]]@coordinates[,c(5,4)]))
Spatial_D00 <- lapply(Spatial_D00, function(x) {colnames(x) <- paste0("Spatial_", 1:2)
                                               x
                                               })  
# Spatial_D00 <- lapply(Spatial_D00, setNames, nm = paste0("Spatial_", 1:2))

lapply(Spatial_D00, function(x) {x$Spatial_2 <- -(x$Spatial_2) 
                                x
                                })
#     return(Spatial_D00)
}
Spatial_D00 <- Spatial_D00_fun(D00)

                      
Spatial_D00_all_fun <- function(Spatial_D00){
    
    tmp <- lapply(seq_along(Spatial_D00), function(i) {bind_rows(Spatial_D00[-i])})
    
    tmp <- lapply(tmp, function(x) {x$Spatial_1<- 0
                             x
                             })
    tmp <- lapply(tmp, function(x) {x$Spatial_2<- 0
                             x
                             })
    
    tmp <- lapply(seq_along(Spatial_D00), function(i) {as.matrix(rbind(Spatial_D00[[i]],tmp[[i]]))
        
        
    })
    
}
           
Spatial_D00_all <- Spatial_D00_all_fun(Spatial_D00) 
                      
temp_fun <- function(D00){

temp <- lapply(D00,function(x) as.data.frame(x@assays[[1]]@counts))
temp <- lapply(temp, function(x)  {x$region <- rownames(x)
                                  x
                                  })  
         lapply(temp, function(x) {rownames(x) <- NULL
                                  x
                                  })      

    }
                
temp <- temp_fun(D00)                      
                      
                          
# merge seurat objects
combined_mat <-reduce(temp, full_join, by = "region");

rownames(combined_mat) <- combined_mat$region
combined_mat$region<- NULL
# removd extra cells
extra_cells <- setdiff(colnames(combined_mat),rownames(Spatial_D00_all[[1]]))
combined_mat <- combined_mat[,which(!colnames(combined_mat)%in%extra_cells)]
combined_mat <- as.matrix(combined_mat)
    
# clean columns of meta data per sample that attached sample's name before rbind
l <- D00
l <- lapply(l, function(x) { colnames(x@meta.data) <- gsub(paste0("_",Assays(x)),"",colnames(x@meta.data));x})
D00 <- l


# first get the list of meta data
list_of_metadata <- lapply(D00, function(x) x@meta.data)
# rbind meta data per samples
meta.data <- do.call("rbind", list_of_metadata)
write.csv(meta.data,'req_meta_data.csv', row.names = T)

    
combined <- CreateSeuratObject(
  counts = combined_mat,
  assay = "scATAC",
  meta.data = meta.data
)

combined@meta.data$Clusters <- factor(combined@meta.data$Clusters,levels = c(
paste0("C",seq_along(unique(combined@meta.data$Clusters)))))

Spatial_D00 <- list()
for (i in seq_along(samples)) {
    Spatial_D00[[i]] <- Spatial_D00_all[[i]][colnames(combined),]
    combined[[paste0(samples[i],"Spatial")]] <- CreateDimReducObject(embeddings = Spatial_D00[[i]]
                                              , key = paste0(samples[i],"Spatial_")
                                              , assay = DefaultAssay(combined))
    
}                           
                           
# we need to run Variable Features
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
                           
UMAPHarmony <- UMAPHarmony[match(colnames(combined),UMAPHarmony$X),]
rownames(UMAPHarmony) <- UMAPHarmony$X
UMAPHarmony$X <- NULL
UMAPHarmony <- as.matrix(UMAPHarmony)

combined[["UMAPHarmony"]] <- CreateDimReducObject(embeddings = UMAPHarmony
                                              , key = "UMAPHarmony_"
                                              , assay = DefaultAssay(combined))
                           
                           
return(combined)
    
    

}


#===========================#===========================#===========================#=================

spatial_in_tissue.obj <- inputs1
#spatial_in_tissue_motif.obj <- readRDS(inputs2)

# creating motif object
find_samples_name <- function(seurat_lst){
  
  sapply(seq_along(seurat_lst), function(i) unique(seurat_lst[[i]]@meta.data$Sample))
  
  
}

samples <- find_samples_name(spatial_in_tissue.obj)
meta.data <- as.data.frame(getCellColData(proj3))
metaData <- list()
for (i in seq_along(samples)) {
  metaData[[i]] <- as.data.frame(getCellColData(proj3[which(proj3$Sample == samples[i])]))
}

motif_object <- list()

for (i in seq_along(samples)) {
  
  motif_object[[i]] <- CreateSeuratObject(counts = dev_score2
                                          , assay = "Spatial", meta.data = metaData[[i]])
  
  DefaultAssay(motif_object = image) <- "Spetial"
  motif_object[[i]][[paste0(samples[i],"_spatial")]] <- spatial_in_tissue[[i]]@images[[1]]
  
}

spatial_in_tissue_motif.obj <- motif_object


UMAPHarmony <- read.csv(inputs3)

combined <- main_func(spatial_in_tissue.obj)
combined_m <- main_func(spatial_in_tissue_motif.obj)


#===========================
scConf1 = createConfig(combined)
makeShinyFiles(combined, scConf1, gex.assay = "scATAC", gex.slot = "data",
               gene.mapping = TRUE, shiny.prefix = "sc1",
          #     shiny.dir = "../../../extra_analysis/ShinyApp/shinyAppMulti_Pieper_Lab/",
               default.gene1 = "Tiam1", default.gene2 = "Ccbe1",
          #       default.multigene = c("DCT","PDCD1","CD19","TIGIT","CD8A","SEMA4D","GZMA","CCL3","CD4","SDF4","B3GALT6"),
               default.dimred = c("UMAPHarmony_1", "UMAPHarmony_2"))




scConf2 = createConfig(combined_m)
makeShinyFiles(combined_m, scConf2, gex.assay = "scATAC", gex.slot = "counts",
              gene.mapping = TRUE, shiny.prefix = "sc2",
 # shiny.dir = "../../../extra_analysis/ShinyApp/shinyAppMulti_Pieper_Lab/",
              default.gene1 = "RFX3-1018", default.gene2 = "NEUROG2-1580",
  # default.multigene = c('WT1-266','SP6-275','SP5-279','SP4-180','TFAP2C-3','ZFX-158','SP2-232','SP1-267',),
              default.dimred = c("UMAPHarmony_1", "UMAPHarmony_2"))


citation = list(
  title   = paste0(project," Data Analysis")
                       )
makeShinyCodesMulti(
  shiny.title = paste0(project,"_Lab Data Analysis"), 
    shiny.footnotes = citation,
  shiny.prefix = c("sc1"
                , "sc2"
                  ),

  shiny.headers = c("Gene Accessibility"
                      , "Peak/Motifs"
                   ), 
  shiny.dir = "./shinyApp"
) 

#==============================================================

sc1def <- readRDS("/root/results/shinyApp/sc1def.rds")
sc2def <- readRDS("/root/results/shinyApp/sc2def.rds")


find_samples_name <- function(lst){

    sapply(seq_along(lst), function(i) unique(lst[[i]]@meta.data$Sample))
}   
samples <- find_samples_name(spatial_in_tissue.obj)


D00 <- list()
for (i in seq_along(samples)) {
D00[[i]] <- spatial_in_tissue.obj[[i]]
nal_cols <- which(colSums(is.na(D00[[i]]@assays[[1]]@counts))>0)
toRemove <- names(nal_cols)
D00[[i]] <- D00[[i]][,!colnames(D00[[i]]) %in% toRemove]
}
Spatial_D00 <- list()
for (i in seq_along(samples)) {
Spatial_D00[[i]] <- as.data.frame(D00[[i]]@images[[1]]@coordinates[,c(5,4)])
colnames(Spatial_D00[[i]]) <- paste0("Spatial_", 1:2)
Spatial_D00[[i]]$Spatial_2 <- -(Spatial_D00[[i]]$Spatial_2)
}
l <- Spatial_D00
xlim <- lapply(l, function(x) { xlim <- c(min(x[,1]),max(x[,1]));xlim})
ylim <- lapply(l, function(y) { ylim <- c(min(y[,2]),max(y[,2]));ylim})




sc1def$limits <- list()
for (i in seq_along(samples)) {
sc1def[['limits']][[paste0(samples[i],"Spatial1")]]<- c(
                                            min(xlim[[i]]),max(xlim[[i]])
                                           ,min(ylim[[i]]),max(ylim[[i]])
                                          )


}




sc2def$limits <- list()
for (i in seq_along(samples)) {
sc2def[['limits']][[paste0(samples[i],"Spatial1")]]<- c(
                                           min(xlim[[i]]),max(xlim[[i]])
                                          ,min(ylim[[i]]),max(ylim[[i]])
                                         )


}




#==================================



sc1def$meta1<- "Clusters"
sc1def$meta2<- "Condition"
sc1def$meta3<- "Sample"

sc2def$meta1<- "Clusters"
sc2def$meta2<- "Condition"
sc2def$meta3<- "Sample"


sc1def$Condition1 <- unique(combined@meta.data$Condition)[1]
sc1def$Condition2 <- unique(combined@meta.data$Condition)[2]

sc2def$Condition1 <- unique(combined_m@meta.data$Condition)[1]
sc2def$Condition2 <- unique(combined_m@meta.data$Condition)[2]



sc1def$Clusters <- read.csv(req_genes1)$x
sc1def$Condition <- read.csv(req_genes2)$x
sc1def$Sample <- read.csv(req_genes3)$x

sc2def$Clusters <- read.csv(req_motifs1)$x
#sc2def$Condition <- read.csv(req_motifs2)$x
#sc2def$Sample <- read.csv(req_motifs3)$x


sc1def$dimred[3] <- paste0(names(combined@reductions)[1],"1")
sc1def$dimred[4] <- paste0(names(combined@reductions)[1],"2")
sc1def$dimred[5] <- paste0(names(combined@reductions)[2],"1")
sc1def$dimred[6] <- paste0(names(combined@reductions)[2],"2")


sc2def$dimred[3] <- paste0(names(combined_m@reductions)[1],"1")
sc2def$dimred[4] <- paste0(names(combined_m@reductions)[1],"2")
sc2def$dimred[5] <- paste0(names(combined_m@reductions)[2],"1")
sc2def$dimred[6] <- paste0(names(combined_m@reductions)[2],"2")


saveRDS(sc1def,"/root/results/shinyApp/sc1def.rds")
saveRDS(sc2def,"/root/results/shinyApp/sc2def.rds")


sc1conf <- readRDS("/root/results/shinyApp/sc1conf.rds")
sc2conf <- readRDS("/root/results/shinyApp/sc2conf.rds")

fav <- which(sc1conf$ID == "Clusters" | sc1conf$ID == "Condition" | sc1conf$ID == "Sample")
rest <- which(!sc1conf$ID == "Clusters" & !sc1conf$ID == "Condition" & !sc1conf$ID == "Sample")

sc1conf <- sc1conf[c(fav,rest),]



fav <- which(sc2conf$ID == "Clusters" | sc2conf$ID == "Condition" | sc2conf$ID == "Sample")
rest <- which(!sc2conf$ID == "Clusters" & !sc2conf$ID == "Condition" & !sc2conf$ID == "Sample")

sc2conf <- sc2conf[c(fav,rest),]


saveRDS(sc1conf,"/root/results/shinyApp/sc1conf.rds")
saveRDS(sc2conf,"/root/results/shinyApp/sc2conf.rds")



rawPath <- "/root/results/"
dataFiles <- dir(rawPath, "*.csv$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)


dataFiles <- dir(rawPath, "*.rds$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

dataFiles <- dir(rawPath, "*.R$", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

dataFiles <- dir(rawPath, "inpMarkers.txt", ignore.case = TRUE, all.files = TRUE)
dataPath <- "/root/results/shinyApp/"
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)


#rawPath <- "/root/"
#dataFiles <- dir(rawPath, "*.R$", ignore.case = TRUE, all.files = TRUE)
#dataPath <- "/root/results/shinyApp/"
#file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)


if (length(unique(proj3$Condition))<=1){
  
  file.remove("/root/results/shinyApp/ui.R")
  file.rename("/root/results/shinyApp/ui_v2.R","/root/results/shinyApp/ui.R")
} else {
  
}
 
 
system ("cp -r /root/results/shinyApp /root/shinyApp")
system ("cp -r /root/results/ArchRProject /root/shinyApp/ArchRProject")

system ("rm -r /root/results")

setwd("/root")
system("ls -a |grep -v 'shinyApp' | xargs rm -r")

file.rename(list.files(pattern="shinyApp"), paste0(project,"_shinyApp"))


