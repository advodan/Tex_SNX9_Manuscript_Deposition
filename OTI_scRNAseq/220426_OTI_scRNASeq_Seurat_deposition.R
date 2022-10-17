# Created by Marcel Trefny in May 2022 for the analysis of OTI T cells scRNASeq of d13 post transfer from tumors of T cells +/- Snx9 KO 
#  Download mm10_ensdb_102_sce.rds from GEO under accession number GSE210535 and place into "input" subfolder 

##############
## Setup
##############

# install.packages('Seurat')
# install.packages('BiocManager')
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("scater")
# BiocManager::install("dittoSeq")
# BiocManager::install('biomaRt')
# BiocManager::install('future')
# BiocManager::install("org.Mm.eg.db", character.only = TRUE)
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("Signac")
# BiocManager::install("SeuratWrappers")
# BiocManager::install("monocle")
# BiocManager::install("Matrix")
# BiocManager::install("patchwork")
# BiocManager::install("dplyr")
# BiocManager::install("ggplot2")
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(dittoSeq)
library(future)
library(Signac)
library(Seurat)
library(monocle)
library(Matrix)
library(ggplot2)
library(patchwork)
library(ggplot2)
library("org.Mm.eg.db", character.only = TRUE)
path = #INSERT PATH HERE

#use paralellization 
plan("multisession", workers = 4)

## Preparation 
 # extracting cell cycle gene homologues from mouse instead of human
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  # Now we extract symbol (external_gene_name), description, gene_biotype (eg whether it is a protein coding or other type of gene):
  ensembl_human_to_mouse <- getBM(attributes=c("external_gene_name", "mmusculus_homolog_associated_gene_name"),
                                  filter="external_gene_name",
                                  values=cc.genes$s.genes,
                                  mart=ensembl)
  # check the structure of the resulting data frame:
  murine_cc_s <- unique(ensembl_human_to_mouse$mmusculus_homolog_associated_gene_name)
  murine_cc_s <- murine_cc_s[nzchar(murine_cc_s) & !is.na(murine_cc_s)]
  write.table(data.frame("symbol" = murine_cc_s), file = "./input/murine_cc_s.txt", sep = "\t", quote = FALSE, row.names = FALSE)

  ensembl_human_to_mouse2 <- getBM(attributes=c("external_gene_name", "mmusculus_homolog_associated_gene_name"),
                                filter="external_gene_name",
                                values=cc.genes$g2m.genes,
                                mart=ensembl)
  # check the structure of the resulting data frame:
  murine_cc_g2m <- unique(ensembl_human_to_mouse2$mmusculus_homolog_associated_gene_name)
  murine_cc_g2m <- murine_cc_g2m[nzchar(murine_cc_g2m)]
  write.table(data.frame("symbol" =murine_cc_g2m), file = "./input/murine_cc_g2m.txt", sep = "\t", quote = FALSE, row.names = FALSE)

################################################################
## Setup
################################################################
setwd(path)
  # check the current active plan
  plan("multisession", workers = 4)
cols = c("#A9A9A9",  "#108080")

# Load the PBMC dataset
  murine_cc_s <- read.table(file = "./input/murine_cc_s.txt",sep = "\t",header = TRUE )$symbol
  murine_cc_g2m <- read.table(file = "./input/murine_cc_g2m.txt",sep = "\t",header = TRUE )$symbol
 
#! #! #! #!#! #!#! #!#! #!#! #!
# mm10_ensdb_102_sce.rds Download this file from GEO under accession number GSE210535 and place in folder
#! #! #! #!#! #!#! #!#! #!#! #!
joined.data <- readRDS("./input/mm10_ensdb_102_sce.rds")

rownames(joined.data) <- rowData(joined.data)$uniqueSymbol
rowData(joined.data)
colData(joined.data)
colData(joined.data)$SampleGroup <- gsub("MT_OTI_Int", "Intergenic", colData(joined.data)$SampleGroup)
colData(joined.data)$SampleGroup <- gsub("MT_OTI_SNX9", "Snx9_KO", colData(joined.data)$SampleGroup)
table(colData(joined.data)$SampleGroup)

grep("ENSM",rownames(cells), value = TRUE)
rownames(joined.data) <- gsub("_ENSM.*", "",rownames(joined.data))

##############
## transform sce to Seurat

# Initialize the Seurat object with the raw (non-normalized data).
cells <- CreateSeuratObject(counts = counts(joined.data),  project = "OTI", min.cells = 3, min.features = 200)
#add info on cells from sce object
cells <- AddMetaData(cells, data.frame(colData(joined.data)), col.name = NULL)
rm(joined.data) #dont need joined data anymore

rownames(cells)
#there are some genes with gene ids in their gene name, clean these

Idents(object = cells) <- cells$SampleGroup
cells[[]]


##############
## QC and Filtering
pdf(file = "./output/plots/QC_filtering.pdf")
  cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^mt-")
  # ^r[S|I|L]
  VlnPlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
  
  VlnPlot(cells, features = c("nFeature_RNA"), ncol = 1) +
    geom_hline(yintercept =  c(2000))
  VlnPlot(cells, features = c("nCount_RNA"), ncol = 1) +
    geom_hline(yintercept =  5000)
  
  plot1 <- FeatureScatter(cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  cells <- subset(cells, subset = nFeature_RNA > 2000  & nCount_RNA>5000 )
  VlnPlot(cells, features = c("percent.mt"), ncol = 1, log = TRUE) +
    geom_hline(yintercept =  15)
  cells <- subset(cells, subset = percent.mt < 15)
  
  #normalize the data
  cells <- NormalizeData(cells, normalization.method = "LogNormalize", scale.factor = 10000)
  
  ## Varible Genes 
  # Identify the 10 most highly variable genes
  cells <- FindVariableFeatures(cells, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(cells), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(cells)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot2
dev.off()

##############
## Scaling and PCA

pdf(file = "./output/plots/PCA_beforeIntegration.pdf")

  cells <- CellCycleScoring(cells, s.features = murine_cc_s, g2m.features = murine_cc_g2m, set.ident = TRUE)
  cells <- AddModuleScore(cells, features = list(grep("^H[1-4]|^Hist", rownames(cells), value = TRUE)) , name = "Histone")
  cells <- AddModuleScore(cells, features = list(grep("^Ifi|^Isg", rownames(cells), value = TRUE)) , name = "Ifi_Isg")
  # 
  cells <- AddModuleScore(cells, features = list(grep("^Rpl|Rps", rownames(cells), value = TRUE)) , name = "Ribosome")
  # 
  all.genes <- rownames(cells)
  cells <- ScaleData(cells, features = all.genes)
  cells <- RunPCA(cells, features = VariableFeatures(object = cells))
  Idents(object = cells) <- cells$Phase
  DimPlot(cells, reduction = "pca")
  FeaturePlot(cells, reduction = "pca", feature = "percent.mt")
  FeaturePlot(cells, reduction = "pca", feature = "nCount_RNA")
  FeaturePlot(cells, reduction = "pca", feature = "nFeature_RNA")
  FeaturePlot(cells, reduction = "pca", feature = "nFeature_RNA", dim = c(2,3))
  FeaturePlot(cells, reduction = "pca", feature = "Histone1")
  FeaturePlot(cells, reduction = "pca", feature = "Histone1", dims = c(2,3))
  FeaturePlot(cells, reduction = "pca", feature = "Ifi_Isg1")
  FeaturePlot(cells, reduction = "pca", feature = "Ifi_Isg1", dims = c(2,3))
  FeaturePlot(cells, reduction = "pca", feature = "Ifi_Isg1", dims = c(3,4))
  # FeaturePlot(cells, reduction = "pca", feature = "Ribosome1", dims = c(1,2))

dev.off()
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# TAKES LONG TO RUN! Therefore skip and start at readRDS below if run before
# with Cell cycle regression, did not include Ribosome, because no effects
cells <- ScaleData(cells, features = all.genes, vars.to.regress = c("S.Score", "G2M.Score", "nFeature_RNA", "Histone1", "Ifi_Isg1"))
cells <- RunPCA(cells, features = VariableFeatures(object = cells))
saveRDS(cells, file = "./rds/cells_integrated_cc_nFeatures.rds")


########################################################################
## Setup from preparated data

setwd(path)

plan("multisession", workers = 4) # use parallelization
cells <- readRDS(file = "./rds/cells_integrated_cc_nFeatures.rds") #load data directly already integrated

pdf(file = "./output/plots/PCA_afterIntegration.pdf")
  
  #Continue with PCA
  print(cells[["pca"]], dims = 1:5, nfeatures = 5)
  VizDimLoadings(cells, dims = 1:2, reduction = "pca")
  VizDimLoadings(cells, dims = 2:3, reduction = "pca")
  Idents(object = cells) <- cells$Phase
  DimPlot(cells, reduction = "pca")
  DimHeatmap(cells, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(cells, dims = 1:15, cells = 500, balanced = TRUE)
  FeaturePlot(cells, reduction = "pca", feature = "percent.mt")
  FeaturePlot(cells, reduction = "pca", feature = "nCount_RNA")
  FeaturePlot(cells, reduction = "pca", feature = "nFeature_RNA")
  
  FeaturePlot(cells, reduction = "pca", feature = "nFeature_RNA", dim = c(2,3))
  FeaturePlot(cells, reduction = "pca", feature = "Histone1")
  FeaturePlot(cells, reduction = "pca", feature = "Histone1", dims = c(2,3))
  FeaturePlot(cells, reduction = "pca", feature = "Ifi_Isg1")
  FeaturePlot(cells, reduction = "pca", feature = "Ifi_Isg1", dims = c(2,3))
  FeaturePlot(cells, reduction = "pca", feature = "Ifi_Isg1", dims = c(3,4))
  
  ElbowPlot(cells, ndims = 30)
dev.off()
#Large number of PCs necessary to describe dataset. Use 20 

##############
## UMAP and Clustering

cells <- FindNeighbors(cells, dims = 1:20)

# UMAP
cells <- RunUMAP(cells, dims = 1:20)

pdf(file = "./output/plots/UMAPs.pdf")
  FeaturePlot(cells, reduction = "umap", feature = "percent.mt")
  FeaturePlot(cells, reduction = "umap", feature = "nCount_RNA")
  FeaturePlot(cells, reduction = "umap", feature = "nFeature_RNA")
  FeaturePlot(cells, reduction = "umap", feature = "Histone1")
  FeaturePlot(cells, reduction = "umap", feature = "Ifi_Isg1")
dev.off()

pdf(file = "./output/plots/UMAPs_clusters_resolution_optim.pdf")
  res <- c(0.3,0.4,0.5,0.6,0.7)
  for(i in 1:length(res)){
    reso <- res[i]
    cells <- FindClusters(cells, resolution = reso)
    print(DimPlot(cells, reduction = "umap", label = TRUE) + ggtitle(paste("Resolution", reso)))
  }
dev.off()

pdf(file = "./output/plots/UMAPs_clusters.pdf")

  #set the resolution here to 0.5 as this appears to generate the most reasonable number of clusters
  cells <- FindClusters(cells, resolution = 0.5)
  Idents(object = cells) <- cells$seurat_clusters
  DimPlot(cells, reduction = "umap", label = TRUE)
  DimPlot(cells, reduction = "umap", split.by = "ExternalSampleName")
  DimPlot(cells, reduction = "umap", group.by = "ExternalSampleName")
  
  DimPlot(cells, reduction = "umap", group.by = "SampleGroup")
  DimPlot(cells, reduction = "umap", group.by = "Phase")
  DimPlot(cells, reduction = "umap", split.by = "Phase", label = TRUE)
  
  FeaturePlot(cells, reduction = "umap", feature = "Ifi_Isg1")
  
  dittoBarPlot( object = cells, var = "Phase",group.by = "SampleGroup")
  dittoBarPlot( object = cells, var = "Phase",group.by = "ExternalSampleName")
  
  DimPlot(cells, reduction = "umap", split.by = "SampleGroup")
  
  genes_of_interest <- c("Ccl5", "Xcl1", "Mki67", "Prf1", "Ccl4", "Nr4a2", "Nr4a3", "Snx9", "Pdcd1", "Havcr2", "Lag3", "Ccr7", "Sell", "Entpd1", "Tcf7", "Tox")
  for(i in 1:length(genes_of_interest)){
    p <- FeaturePlot(cells, reduction = "umap", feature = genes_of_interest[i])
    print(p)
  }
  
  FeaturePlot(cells, reduction = "umap", feature = "Satb1")
  FeaturePlot(cells, reduction = "umap", feature = "Irf8")
  FeaturePlot(cells, reduction = "umap", feature = "Runx3")
  FeaturePlot(cells, reduction = "umap", feature = "Tgfbr1")
  FeaturePlot(cells, reduction = "umap", feature = "Il2ra")
  FeaturePlot(cells, reduction = "umap", feature = "Il2rb")
  FeaturePlot(cells, reduction = "umap", feature = "Ncam1")

dev.off()


########## 
## Signatures
#load signatures
signatures <- sapply(list.files("./signatures/", pattern = ".txt"), function(x) read.table(paste0("./signatures/", x), header = FALSE, sep = "\t")[,1], simplify = FALSE, USE.NAMES = TRUE )

#clean names
names(signatures) <- gsub(".txt", "", names(signatures))
signatures

#add score to each cell from each signature
# names(signatures) <- sapply(names(signatures), function(x) {gsub(paste0(x,"1"), x, names(signatures)); x} )

cells <- AddModuleScore(cells, features = signatures, name = names(signatures))

#the AddModuleScore has an amazingly stupid feature that adds a number to each column. Thus loop over all the signatures and remove the numbering...
# run three times (space for 30 signatures)
for(k in 1:3){
for(i in 1:length(signatures)){
  sign <- names(signatures)[i]
  colnames(cells@meta.data) <-
    gsub(x = colnames(cells@meta.data)
         , pattern = paste0(sign,"[1-9].*")
         , replacement = sign)
}
}
## plot each signature onto umap and a violin plot
#function to generate mean +/- 1sd in violin plot
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


pdf(file = "./output/plots/Signatures_dotplots.pdf", width = 5, height = 5)
    # a Dotplot for each group of signatures for all clusters and signatures
    Idents(object = cells) <- cells$seurat_clusters
    signature_groups <- c("Andreatta",  "Schietinger", "Miller")
    for(i in 1:length(signature_groups)){
        p <- DotPlot(cells, feature = names(signatures)[grepl(signature_groups[i], names(signatures) )] ) + scale_color_gradient2(low = "blue", high = "red", mid = "white")  +
        theme(axis.text.x = element_text(angle = -45, hjust = 0))
      print(p)
    }
dev.off()

pdf(file = "./output/plots/Signatures.pdf", width = 5, height = 5)
    for(i in 1:length(signatures)){
      sign <- names(signatures)[i]
      p <- FeaturePlot(cells, reduction = "umap", feature = sign) 
      print(p)
      
    }
    # a violin comparing the different clusters
    Idents(object = cells) <- cells$seurat_clusters
    for(i in 1:length(signatures)){
      sign <- names(signatures)[i]
      p <- VlnPlot(cells, features = sign, pt.size = 0) +  stat_summary(fun.data="data_summary", 
                                                                        geom="pointrange", color="black")
      print(p)
    }
    
    #plot cell cycle per cluster
    p1 <- VlnPlot(cells, features = "S.Score", pt.size = 0) +  stat_summary(fun.data="data_summary", 
                                                                      geom="pointrange", color="black")
    p2 <-  VlnPlot(cells, features = "G2M.Score", pt.size = 0) +  stat_summary(fun.data="data_summary", 
                                                                             geom="pointrange", color="black")
    print(p1+p2)
    
    Idents(object = cells) <- cells$SampleGroup
    #now a violin comparing integenic vs KO
    for(i in 1:length(signatures)){
      sign <- names(signatures)[i]
      p <- VlnPlot(cells, features = sign, pt.size = 0) +  stat_summary(fun.data="data_summary", 
                                                                        geom="pointrange", color="black")
      print(p)
    }
    #plot cell cycle per cell type
    p1 <- VlnPlot(cells, features = "S.Score", pt.size = 0) +  stat_summary(fun.data="data_summary", 
                                                                            geom="pointrange", color="black")
    p2 <-  VlnPlot(cells, features = "G2M.Score", pt.size = 0) +  stat_summary(fun.data="data_summary", 
                                                                               geom="pointrange", color="black")
    print(p1+p2)

dev.off()

###################################
# rename clusters and based on these names plot umap for figures

Idents(object = cells) <- cells$seurat_clusters
new.cluster.ids <- c("Tex-term", "Tex-prolif", "Tex-Prf1", "Tem", "Tprolif", "Tpex",
                     "Tex-early-prolif", "Tnaive-like", "Tcytokine")
names(new.cluster.ids) <- levels(cells)
cells <- RenameIdents(cells, new.cluster.ids)
cells$cluster_names <- Idents(cells)
cells$cluster_names <- factor(cells$cluster_names, levels = new.cluster.ids )
cells[[]]
cells$umap1 <- cells@reductions$umap@cell.embeddings[,"UMAP_1"]
cells$umap2 <- cells@reductions$umap@cell.embeddings[,"UMAP_2"]

saveRDS(cells, file = "./rds/cells_integrated_clustered_annotated.rds")

#### late entry point
setwd(path)
# check the current active plan
plan("multisession", workers = 4)
cols = c("#A9A9A9",  "#108080")
cells <- readRDS(file = "./rds/cells_integrated_clustered_annotated.rds")

cells[[]] #metadata
write.table(cells[[]], file = "output/lists/OTI_scRNAseq_metadata_cells.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

pdf(file = "./output/plots/UMAPs_clusters_Figures.pdf", width = 7, height = 4)
  Idents(object = cells) <- cells$cluster_names
  DimPlot(cells, reduction = "umap", split.by = "SampleGroup") 
dev.off()


# a Dotplot for each group of signatures for all clusters and signatures
pdf(file = "./output/plots/Signatures_dotplots_renamed.pdf", width = 5.1, height = 5.5)
  Idents(object = cells) <- cells$cluster_names
  signature_groups <- c("Andreatta",  "Schietinger", "Miller")
  for(i in 1:length(signature_groups)){
    p <- DotPlot(cells, feature = names(signatures)[grepl(signature_groups[i], names(signatures) )], cols = "RdBu" ) +
      theme(axis.text.x = element_text(angle = -45, hjust = 0), legend.text =  element_text(size = 10), legend.title =  element_text(size = 12)) +
      xlab("")
    print(p)
    
    d <- p$data
    write.table(d, file = paste0("output/lists/signatures_dotplots_", signature_groups[i],".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
dev.off()

##############
## Markers

Idents(object = cells) <- cells$cluster_names
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top.markers <- cells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>% dplyr::select(gene)
top.markers <- top.markers$gene

write.table( cells.markers, file = "./output/lists/markers_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)


pdf(file = "./output/plots/marker_plots.pdf", width = 9, height = 5)
  
  ## Markers Clusters
  #plot markers per cluster
  Idents(object = cells) <- cells$cluster_names
  FeaturePlot(cells, reduction = "umap", feature = top.markers)
  #selected exhaustion and function
  p <- DotPlot(cells, feature = c("Tox",  "Lag3", "Havcr2","Nkg7", "Gzmb","Gzmc","Prf1",  "Cxcr6", "Cxcr3","Ly6c2", "Mki67","Ccnb2", "Nme1", "Tcf7","Bach2", "Batf", "Slamf6", "Lef1", "Ccr7"), cols = "RdBu") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),axis.text.y = element_text(size = 14), legend.text =  element_text(size = 10), legend.title =  element_text(size = 12)) + xlab("") 
  print(p)
  d <- p$data
  write.table(d, file = "output/lists/markers_dotplot_1.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  #functionality
  DotPlot(cells, feature = c(  "Il2", "Ifng",  "Gzmb", "Gzmc","Prf1", "Ccl3", "Ccl4", "Ccl5"), cols = "RdBu") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text =  element_text(size = 10), legend.title =  element_text(size = 12)) + xlab("") 
  #transcription factors
  p <- DotPlot(cells, feature = c("Tox", "Batf", "Runx3", "Eomes", "Tbx21", "Nfatc1","Nfatc2","Nr4a1" ,"Nr4a2","Nr4a3","Tcf7", "Irf4", "Irf8", "Id1", "Id2", "Id3", "Bach2", "Bcl6"), cols = "RdBu") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text =  element_text(size = 10), legend.title =  element_text(size = 12)) + xlab("") 
  print(p)
  d <- p$data
  write.table(d, file = "output/lists/markers_dotplot_TF.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  #find markers for each of the clusters
  cells.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5
  p <- DoHeatmap(cells,WhichCells(cells, downsample = 500, seed = 1), features = top5$gene) + NoLegend()
  print(p)
  d <- p$data
  write.table(d, file = "output/lists/markers_heatmap_large_ds500.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  ## KO vs intergenic
  #for each cluster show the proportion of cells per conditions
  #this plot is interesting but it doesnt account for differences in cell numbers per sample. thus use the proportions below
  dittoBarPlot( object = cells, var = "cluster_names",group.by = "SampleGroup")
  dittoBarPlot( object = cells, var = "SampleGroup",group.by = "cluster_names") + scale_fill_manual(values = cols)
  
dev.off() 
  
## Proportion of cells in each cluster
pdf(file = "./output/plots/proportions_plots.pdf", width = 8, height = 7)
  
  tab <- table(data.frame(cells$SampleGroup, cells$cluster_names))
  tab.rowsums <- rowSums(tab)
  proportions_clusters <- data.frame(t(tab/tab.rowsums))
  
  proportions_clusters
  
  ggplot(proportions_clusters, aes(x = cells.cluster_names, y = Freq, fill = cells.SampleGroup)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle("Cluster Proportions") +
    scale_fill_manual(values = cols, aesthetics = c("colour", "fill")) +
    xlab("") +
    ylab("Frequency of cells among condition") +
    labs(fill = "Condition") +
    theme_classic() +
    theme(text = element_text(size = 15), axis.text =   element_text(size = 15), axis.text.x =   element_text(angle = 45, hjust= 1))
  #write to file
write.table(proportions_clusters, file = "./output/lists/proportions_clusters.txt", sep = "\t", row.names = FALSE, quote = FALSE )
  
  Idents(object = cells) <- cells$SampleGroup
  VlnPlot(cells, features = c("Ccl5", "Nr4a2", "Nkg7", "Ccl4", "Ccl3", "Il10","Cxcr3", "Irf8", "Fosb"), pt.size = 0)
  
  Idents(object = cells) <- cells$cluster_names
  DotPlot(subset(x= cells, idents = "Tex-term"), features = c("Nr4a2", "Tox" ,"Ccl5", "Ccl4", "Ccl3", "Cxcr3"), group.by =  "SampleGroup")+
    coord_flip()
    
dev.off()

# marker plots as tiff
tiff(file = "./output/plots/Markers_heatmap_unnamed.tiff", width = 1000, height = 900, res = 170)
  Idents(object = cells) <- cells$seurat_clusters
  DoHeatmap(cells,WhichCells(cells, downsample = 500, seed = 1), features = top5$gene) + NoLegend()
dev.off()

tiff(file = "./output/plots/Markers_heatmap_named.tiff", width = 1000, height = 900, res = 170)
  Idents(object = cells) <- cells$cluster_names
  DoHeatmap(cells,WhichCells(cells, downsample = 500, seed = 1), features = top5$gene) + NoLegend()
dev.off()

# proportion of cells in each cluster
pdf(file = "./output/plots/proportions_figures.pdf", width = 6, height = 5)
  
  ggplot(proportions_clusters, aes(x = cells.cluster_names, y = Freq, fill = cells.SampleGroup)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle("Cluster Proportions") +
    scale_fill_manual(values = cols, aesthetics = c("colour", "fill")) +
    xlab("") +
    ylab("Frequency of cells") +
    labs(fill = "Condition") +
    theme_classic() +
    theme(text = element_text(size = 15), axis.text =   element_text(size = 15, color = "black"), axis.text.x =   element_text(angle = 45, hjust= 1))
dev.off()

#####################################################
## Differential Gene Expression

# dif genes between KO and intergenic in all cells (disregarding the clusters)
pdf(file = "./output/plots/diffGenes_AllCells.pdf")
  Idents(object = cells) <- cells$SampleGroup
  Snx9_diffGenes_allClusters <- FindMarkers(cells, ident.1 = "Snx9_KO", ident.2 = "Intergenic", min.pct = 0.3, logfc.threshold = -1)
  head(Snx9_diffGenes_allClusters, n = 100)
  Snx9_diffGenes_allClusters$gene <- rownames(Snx9_diffGenes_allClusters)
  p <- ggplot(Snx9_diffGenes_allClusters, aes(x= avg_log2FC, y= -log10(p_val_adj), label = gene)) + 
    geom_point()+ 
    geom_label_repel(data = subset(Snx9_diffGenes_allClusters,abs(avg_log2FC) > 0.25 & p_val_adj < 0.01  )) +
    ggtitle("Differentially expressed genes between Snx9 KO and intergenic\n") +
    ylab("-log10(adjusted p-value)") +
    xlab("average log2 fold change")+
    theme_bw()
  print(p)
    #also for differences between KO and integenic make a Heatmap for marker genes
  DoHeatmap(cells, WhichCells(cells, downsample = 50, seed = 1), features = rownames(Snx9_diffGenes_allClusters) ) + NoLegend()
  write.table(Snx9_diffGenes_allClusters, file = paste0("./output/lists/diffGenes_allCells.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

dev.off()

saveRDS(Snx9_diffGenes_allClusters, file = "./output/Snx9_diffGenes_allClusters.rds")

## ### ## ## ### ## ## ### ## ## ### ## ## ### ## 
# TAKES LONG TO RUN! Start below at readRDS for quick start if performed already
# between Snx9 KO and intergenic for each cluster individually
# calculate all diff genes first
Snx9_diffGenes_clusters_aggregated <- data.frame()
Idents(object = cells) <- cells$SampleGroup
  for(j in 1:length(unique(cells$cluster_names))){
        cluster = as.character(unique(cells$cluster_names))[j]
      Snx9_diffGenes_cluster <- FindMarkers(subset(x = cells, subset = cluster_names == cluster ), 
                                            ident.1 = "Snx9_KO", ident.2 = "Intergenic", min.pct = 0.3,
                                            logfc.threshold = -1)
      Snx9_diffGenes_cluster$gene <- rownames(Snx9_diffGenes_cluster)
      Snx9_diffGenes_cluster$cluster <- cluster
      Snx9_diffGenes_clusters_aggregated <- rbind(Snx9_diffGenes_clusters_aggregated,Snx9_diffGenes_cluster )
      write.table(Snx9_diffGenes_cluster, file = paste0("./output/lists/diffGenes_Cluster_", cluster, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
#combine with diff genes of all cells
Snx9_diffGenes_clusters_aggregated <- rbind(Snx9_diffGenes_clusters_aggregated, data.frame(Snx9_diffGenes_allClusters, cluster = "all cells"))
write.table(Snx9_diffGenes_clusters_aggregated, file = paste0("./output/lists/Snx9_diffGenes_clusters_aggregated.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(Snx9_diffGenes_clusters_aggregated, file = "./output/Snx9_diffGenes_clusters_aggregated.rds")
Snx9_diffGenes_clusters_aggregated <- readRDS(Snx9_diffGenes_clusters_aggregated, file ="./output/Snx9_diffGenes_clusters_aggregated.rds" )
head(Snx9_diffGenes_clusters_aggregated, n = 100)

# QUICK ENTRY POINT -> takes long to perform otherwise below
Snx9_diffGenes_clusters_aggregated <- readRDS( file = "./output/Snx9_diffGenes_clusters_aggregated.rds")
head(Snx9_diffGenes_clusters_aggregated, n = 100)

#now plot these differences
pdf(file = "./output/plots/diffGenes_eachCluster.pdf", width = 6, height = 5)
for(j in 1:length(unique(Snx9_diffGenes_clusters_aggregated$cluster))){
  clu = as.character(unique(Snx9_diffGenes_clusters_aggregated$cluster))[j]
  diffGenes <- subset(Snx9_diffGenes_clusters_aggregated, cluster == clu)
  p <- ggplot(diffGenes, aes(x= avg_log2FC, y= -log10(p_val_adj), label = gene, col = avg_log2FC,  size = abs(avg_log2FC))) +
    geom_point()+
    scale_size(range = c(0.1,3)) +
    scale_color_viridis_c(option = "viridis") +
    geom_text_repel(data = subset(diffGenes,abs(avg_log2FC)  > 0.3 & -log10(p_val_adj) > 2 |  -log10(p_val_adj) > 30), col = "black", size = 4.5) +
    ggtitle(paste0("Differential gene expression\n", clu)) +
    ylab("-log10 (adj. p-value)") +
    xlab("log2 fold change Snx9 KO vs intergenic")+
    labs(col = "log2 fold change", size = "abs(log2 fold change)") +
    theme_bw() +
    theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
  print(p)
  }
  
dev.off()

#########
# Selected Differential genes for figures shown as dotlots 
pdf(file = "./output/plots/diffGenes_eachCluster_dotplot.pdf", height = 4, width = 5.2)
  genes_of_interestSign <- c("Snx9","Il2ra","Ccl3", "Ccl4", "Ccl5","Xcl1",  "Gzmb","Cxcr6",   "Nr4a2", "Nr4a1", "Nr4a3", "Ccr7",  "Tox" )
  
  dat <- Snx9_diffGenes_clusters_aggregated %>% filter(gene %in% genes_of_interestSign & p_val_adj < 0.05  ) 
  dat$gene <- factor(dat$gene, levels = genes_of_interestSign)  
  
  write.table(dat, file = "./output/lists/source_data_dotplot_diffGenes.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  ggplot(dat, aes(x = gene, y = cluster, size = -log10(p_val_adj), col = avg_log2FC )) +
    geom_point() +
    theme_bw() +
    scale_size(range = c(2,6) ) +
    labs(col = "log2 fold change", size = "-log10(adj. p value)") +
    scale_color_distiller(palette = "RdBu", limits = c(-1, 1)) +
 #    scale_color_gradient2(low = "blue", high = "red", mid = "white", breaks = c(1,-1)) +
    ggtitle("Selected differentially expressed genes")+
    theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black") , axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.text = element_text(size = 13), legend.title = element_text(size = 13))
dev.off()


