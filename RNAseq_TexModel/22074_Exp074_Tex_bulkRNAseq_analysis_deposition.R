# analysis of the bulk RNAseq data from the Tex model Exp074. This time containing Trest, Ttumor, Teff and Tex from 5 donors at day 12. 

#initialization
# BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
# BiocManager::install('circlize')
# BiocManager::install('edgeR')
# BiocManager::install('ComplexHeatmap')
# BiocManager::install('tidyr')
# BiocManager::install('paletteer')
library(limma)
library(edgeR)
library(locfit)
library(RColorBrewer)
library(grDevices)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(viridis)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(stringr)
library(circlize)
library(paletteer) 

myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
condition_colors <- c("#F8F9F9", "#BFBFC1",  "#FB0207", "#800002")
donor_colors <- c("HD276" = "#89E8ACFF", "HD280" = "#67DBA5FF", "HD286" = "#4CC8A3FF", "HD287" = "#38B2A3FF")
set.seed(123)

#setup environment
setwd(#Insert Path of Script here#)
if(!dir.exists("output/")){dir.create("output")}
if(!dir.exists("output/lists")){dir.create("output/lists")}
if(!dir.exists("output/plots")){dir.create("output/plots")}
if(!dir.exists("output/lists/Zhang")){dir.create("output/lists/Zhang")}
if(!dir.exists("output/lists/Zhang")){dir.create("output/lists/Zhang")}
if(!dir.exists("output/lists/Zhang_Heatmaps")){dir.create("output/lists/Zhang_Heatmaps")}
if(!dir.exists("output/lists/GO_BP")){dir.create("output/lists/GO_BP")}

## Load Count Data and perform some summary statistics and plots

# can also be found under GEO  accession number GSE210534
dgelist <- readRDS("input/ensdb_105_dge_list_final.rds")
dgelist


colnames(dgelist$counts) <- dgelist$samples$ExternalSampleName

dim(dgelist$counts)
colSums(dgelist$counts)
colMeans(dgelist$counts)

hist(log10(rowSums(dgelist$counts)))

#exclude transcripts which are not expressed at least in 1 condition at cpm > 0. 
cpms <- cpm(dgelist$counts, prior.count = 8, log=TRUE)
keep <- (rowSums(cpms>0) >= 4 & !is.na(dgelist$genes$SYMBOL))
summary(keep) 
# Mode   FALSE    TRUE 
# logical   44938   17916 
y <- dgelist[keep,] #only keep these transcripts

#also exclude pseudogenes and small RNA transcripts (mRNA sequencing was performed)
y <- y[(y$genes$METABIOTYPE == "protein_coding" & !is.na(y$genes$METABIOTYPE) ),]
dim(y)
#12398 16
table(y$genes$METABIOTYPE)

#get rid of duplicates of same gene symbol, keep the transcript with the higher expression
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$SYMBOL)
summary(d)
# Mode   FALSE    TRUE 
# logical   11929     469 
y <- y[!d,]
nrow(y) #11929

saveRDS(y, file = "./output/dgeList_Exp074_RNAseq.rds")

#normalize for library sizey
y <- calcNormFactors(y)

#annotate samples
y$samples
y$samples$SampleGroup  <- factor(y$samples$SampleGroup, levels = c("Trest", "Ttumor", "Teff", "Tex"))
rownames(y$counts) <- y$genes$SYMBOL

par(mfrow=c(1,1))
cpms_filtered <- cpm(y, prior.count = 8, log=TRUE, normalized.lib.sizes = T) #log cpms corrected
rpkms_filtered <- rpkm(y, prior.count = 8, log=TRUE, normalized.lib.sizes = T) #log RPKM corrected
rownames(rpkms_filtered) <- y$genes$SYMBOL
rm(cpms)

#compute PCA
summary(apply(cpms_filtered, 1, var) == 0) ## Any gene with no variability?
pca1 <- prcomp(t(cpms_filtered), scale = T)
summary(pca1) ## 59% variance on PCs 1:3

# loadings: coordinates of PCs in original coordinate system
loadings <- pca1$rotation
## scores: coordinates of scaled data (observations) projected onto the PCs.
scores <- pca1$x
plotPCA <- function(n=1, m=2) {
  col.v <- condition_colors
  pchs <- c(21,22)[as.integer(y$samples$sex)]
  pchs <- 21
  pdf(paste0("./output/plots/PCA_", n,"vs", m ,".pdf"), width=7)
  {
    par(mar=c(5,5, 5,5) )
    plot(scores[,n], scores[,m], xlim=c(min(scores[,n])-30, max(scores[,n])+30), 
         ylim=c(min(scores[,m])-30, max(scores[,m])+30),
         xlab=paste("PC", n, ": ", round(summary(pca1)$importance[2,n],3)*100, "% variance explained", sep=""), 
         ylab=paste("PC", m, ": ", round(summary(pca1)$importance[2,m],3)*100, "% variance explained", sep=""),
         col="black", bg=col.v, lwd = 0.5, cex= 2, cex.lab=2, cex.main = 2, 
         pch=pchs, main = "Principal Component Analysis")
    legend("topright", levels(y$samples$SampleGroup), pch=21, cex=1.8, col="black", pt.bg=col.v)
    #legend("bottomright", levels(y$samples$sex), cex=0.95, col="black")}
    dev.off()
  }
}
plotPCA(1,2)
plotPCA(2,3)
plotPCA(1,3)

#save PCA results
write.table(scores, file= "output/lists/PCA_scores.txt", sep= "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(summary(pca1)$importance, file= "output/lists/PCA_importance.txt", sep= "\t", quote = FALSE, row.names = TRUE, col.names = NA)

#for genes of interest make boxplots of expression across conditions
genes_of_interest <- c("PDCD1", grep("SIGLEC", y$genes$SYMBOL, value = TRUE ))
genes_of_interest <- c("PDCD1", "CTLA4", "TOX", "TOX2", "NR4A1", "NR4A2", "NR4A3", "LAG3","PTPN6", "GZMB",   "TNFRSF9", "TCF7")
cpms_genes_of_interest <- data.frame(cpms_filtered[rownames(cpms_filtered) %in% genes_of_interest,])
cpms_genes_of_interest$SYMBOL <- rownames(cpms_genes_of_interest)
cpms_genes_of_interest <- cpms_genes_of_interest  %>% 
                              gather("SampleGroup", "cpm", -"SYMBOL") %>%
                              mutate(donor = gsub("_.*", "", SampleGroup),
                                    condition = gsub(".*_", "", SampleGroup) )
cpms_genes_of_interest$condition <- factor(cpms_genes_of_interest$condition, levels = c("Trest", "Ttumor", "Teff", "Tex"))
pdf(file = "./output/plots/selected_genes_boxplots.pdf", width = 9)
ggplot(cpms_genes_of_interest, aes(x = condition, y = cpm, fill = condition)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center') + facet_wrap(~SYMBOL,  scales = "free") + ylab("log cpm") +
  theme_bw() +
  scale_fill_manual(values = condition_colors)
dev.off()

write.table(cpms_genes_of_interest, file= "output/lists/boxplot_selected_genes_cpms.txt", sep= "\t", quote = FALSE, row.names = TRUE, col.names = NA)


## differential gene expression analysis
#Build Design Matrix, make a batch correction for donor
Donor <- y$samples$Donor
SampleGroup <- y$samples$SampleGroup
design_donor <- model.matrix(~ Donor + SampleGroup)
rownames(design_donor) <- colnames(y)
colnames(design_donor) <- gsub("SampleGroup", "", colnames(design_donor))
design_donor

#Build GLM
y_donor <- estimateDisp(y, design_donor)
plotBCV(y_donor, main="BCV Plot paired") #Plotting Biological Coefficient of Variation

fit_donor <- glmQLFit(y_donor, design_donor, prior.count=8)

#diff gene cutoffs
pcut=0.01
logfcut=0.75

#ANOVA like analysis between the conditions
anova <- glmQLFTest(fit_donor, coef=c(5:ncol(design_donor)) )
anova_all <-  topTags(anova, adjust.method = "BH",  n = Inf, sort.by = "PValue")$table
anova_all


#find only significantly changed genes with a certain logfc change in at least one condition
anova_all.signifGenes <- anova_all[anova_all$FDR < pcut &
                                     (abs(anova_all$logFC.Ttumor) > logfcut | abs(anova_all$logFC.Teff) > logfcut |abs(anova_all$logFC.Tex) > logfcut )    , ]
dim(anova_all.signifGenes) #4242

###################
# Heatmaps
####################
#plot these genes in a heatmap for each sample
pdf(file = "./output/plots/differential_genes_anova_heatmap.pdf", width = 7, height = 6)
  genes_of_interest <- cpms_filtered[rownames(cpms_filtered) %in% anova_all.signifGenes$SYMBOL, ]
  genes_of_interest <- genes_of_interest[,c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)]

  #annotate samples with condition names and colors
  topAn <- HeatmapAnnotation(
    df = data.frame(group = c(rep(1,4),rep(2,4),rep(3,4),rep(4,4) ) ),  
    col = list(group = c("1" = "#F8F9F9", "2"= "#BFBFC1","3"=  "#FB0207", "4"= "#800002") ),
    annotation_legend_param = list(group = list(labels=levels(y_donor$samples$SampleGroup), title="Condition")),
    annotation_name_gp = gpar(fontsize = 0)
  )

  #label only selected genes
  selected_genes <- c("PDCD1", "HAVCR2", "LAG3", "TOX", "TOX2", "NR4A1",  "NR4A3",  "TNFRSF9",  "IFNG", "PTPN6",   "TNF", "RORA",  "TCF7",  "IL7R", "LILRB3")
  rownames(genes_of_interest)[!rownames(genes_of_interest) %in% selected_genes ] <- ""
    nr_clusters <- 6
    set.seed(1234)
    ht <- Heatmap(t(scale(t(genes_of_interest), scale=TRUE, center=TRUE)), col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)), name = "z-score log(cpm)",
                  column_title =  paste("RNAseq Differentially Expressed Genes") , column_title_gp = gpar(fontsize=14), cluster_columns  = FALSE, cluster_rows = TRUE,  
                  row_km = nr_clusters,show_column_names = FALSE,show_row_names = TRUE,show_row_dend = FALSE, show_column_dend = FALSE,
                  row_names_gp = gpar(fontsize = 8.2, justify = "left"),top_annotation = topAn, row_title_gp = gpar(fontsize = 12),
                  heatmap_legend_param = list(labels_gp = gpar(fontsize = 12 ),
                                              title_gp = gpar(fontsize = 12 )), use_raster = FALSE) 
    clusterlist = row_order(ht)
    set.seed(1234)
    draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
    
    #format table for export and save
    genes_of_interest <- cpms_filtered[rownames(cpms_filtered) %in% anova_all.signifGenes$SYMBOL, ]
    genes_of_interest <- genes_of_interest[,c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)]
    expt <-  data.frame(t(scale(t(genes_of_interest), scale=TRUE, center=TRUE)))
    head(expt)      
    expt$SYMBOL <- row.names(expt)
    clusters <- data.frame()
    for(i in 1:length(names(clusterlist))){
      clu <- names(clusterlist[i])
        clusters <- rbind(clusters, data.frame("SYMBOL" = rownames(genes_of_interest[clusterlist[[clu]], ]), "cluster" = clu))
      }
    head(clusters)
    tail(clusters)
    expt <- left_join(expt, clusters, by = "SYMBOL")
    write.table(expt, file = "output/lists/anova_heatmap_scaled_centered_cpms.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    
dev.off()



#write to file all signficant anova genes including their cluster annotation
genes_of_interest <- cpms_filtered[rownames(cpms_filtered) %in% anova_all.signifGenes$SYMBOL, ]
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(SYMBOL = rownames(genes_of_interest[clusterlist[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)
write.table(clu_df, file= "./output/lists/anova_clusters_genelist.txt", sep  = "\t", quote = FALSE, row.names = FALSE)
clu_df[clu_df$SYMBOL == "SNX9",]
clu_df[clu_df$Cluster == "cluster1", ]



############### 
# plots of genes of interest
genes_of_interest <- cpms_filtered[rownames(cpms_filtered) %in% selected_genes, ]
genes_of_interest <- genes_of_interest[selected_genes,] #reorder to cluster conditions
genes_of_interest <- genes_of_interest[,c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)]

topAn <- HeatmapAnnotation(
  df = data.frame(group = c(rep(1,4),rep(2,4),rep(3,4),rep(4,4) ) ),  
  col = list(group = c("1" = "#F8F9F9", "2"= "#BFBFC1","3"=  "#FB0207", "4"= "#800002") ),
  annotation_legend_param = list(group = list(labels=levels(y_donor$samples$SampleGroup), title="Condition")),
  annotation_name_gp = gpar(fontsize = 0)
)

ht <- Heatmap(t(scale(t(genes_of_interest), scale=TRUE, center=TRUE)), col = viridis(256, option = "magma"), name = "z-score log(cpm)",
              column_title =  paste("Selected Genes") , column_title_gp = gpar(fontsize=14), cluster_columns  = FALSE, cluster_rows = FALSE,  
              show_column_names = FALSE,show_row_names = TRUE,show_row_dend = FALSE, 
              row_names_gp = gpar(fontsize = 11, justify = "left"),top_annotation = topAn,
              heatmap_legend_param = list(labels_gp = gpar(fontsize = 12 ),
                                          title_gp = gpar(fontsize = 12 )), use_raster = FALSE)
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))


###################
# Contrasts between conditions
####################
#define contrasts and perform differential gene expression analysis between the conditions
design_donor 
mycontrasts <- list(c(0,0,0,0,1,0,0), c(0,0,0,0,0,1,0), c(0,0,0,0,0,0,1), c(0,0,0,0,-1,1,0), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1) )
names(mycontrasts) <- c("Ttumor_vs_Trest", "Teff_vs_Trest", "Tex_vs_Trest",   "Teff_vs_Ttumor", "Tex_vs_Ttumor", "Tex_vs_Teff")

lrt_all <- sapply(names(mycontrasts),function(x) glmQLFTest(fit_donor,contrast=mycontrasts[[x]] ), simplify=FALSE, USE.NAMES = TRUE )

tables_all_unsorted <- sapply(names(lrt_all), function(x)  topTags(lrt_all[[x]],adjust.method="BH", n = Inf, sort.by = "none")$table,simplify=FALSE, USE.NAMES = TRUE)
tables_all <- sapply(names(lrt_all), function(x)  topTags(lrt_all[[x]],adjust.method="BH",  n = Inf, sort.by = "PValue")$table, simplify=FALSE, USE.NAMES = TRUE)
tables_signif <- sapply(names(lrt_all), function(x) tables_all[[x]][abs(tables_all[[x]]$logFC) > logfcut & tables_all[[x]]$FDR < pcut, ],simplify = FALSE, USE.NAMES = TRUE )


sapply(tables_all, function(x) head(x, n=20),simplify=FALSE, USE.NAMES = TRUE)
sapply(names(lrt_all), function(x) table(decideTestsDGE(lrt_all[[x]], p.value = pcut, adjust.method="BH", lfc = logfcut)), simplify = TRUE, USE.NAMES = TRUE)
de_all <- sapply(names(lrt_all), function(x) decideTestsDGE(lrt_all[[x]], p.value = pcut, adjust.method="BH", lfc = logfcut), simplify = FALSE, USE.NAMES = TRUE)
sapply(names(tables_all), function(x) write.table(tables_all[[x]],paste0("./output/lists/", x,"_dgelist_full.txt"), sep ="\t", quote=FALSE, row.names = FALSE))

#########
# perform enrichment analysis with Zhang Pan Cancer T cell states 
##########
#load gene sets 
customGeneSets_Zhang <- list()
for (f in list.files("./input/custom_genesets/scRNAseq_ZhangPanCancer/", pattern = "*.txt")) {
  
  myGeneSet <- read.table(paste0("./input/custom_genesets/scRNAseq_ZhangPanCancer/", f), h=TRUE, sep="\t")
  cat("  ", length(myGeneSet[,1]), "genes in geneset ")
  myGeneSet <- myGeneSet[!is.na(myGeneSet$ES),] #remove empty lines
  ## If you read the genes as symbols in the files, match them to y$genes abd use the ENTREZID column
  ## Add gene set to list
  clean_names <- gsub(".txt", "", f)
  clean_names <- gsub("CD8_c[0-9][0-9]_", "", clean_names)
  customGeneSets_Zhang[[clean_names]] <- myGeneSet
  cat(clean_names, "\n")
}
length(customGeneSets_Zhang) #17
sapply(customGeneSets_Zhang, function(x) head(x), simplify = FALSE) #show top genes per set

customGeneSets_Zhang <- sapply(customGeneSets_Zhang, function(x) head(x, n = 100), simplify = FALSE) #only use top 100 genes based on effect size ranking (to account for differerent sizes of gene sets)
write.table(names(customGeneSets_Zhang), file = "./output/lists/customGeneSets_Zhang_names.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names =  FALSE)
sapply(names(customGeneSets_Zhang), function(x) write.table(data.frame(customGeneSets_Zhang[[x]]), paste0("./output/lists/Zhang/customgenesets_Zhang_", x, ".txt") ,sep='\t', quote = FALSE ), simplify = FALSE, USE.NAMES = TRUE)



## Custom Gene Sets Enrichment
# create index
index_custom_Zhang <- sapply(names(customGeneSets_Zhang), function(x) y_donor$genes$SYMBOL %in% customGeneSets_Zhang[[x]][,"SYMBOL"] , simplify=FALSE, USE.NAMES=TRUE)
index_custom_Zhang_number <- sapply(names(index_custom_Zhang),function(x) sum(index_custom_Zhang[[x]], na.rm=TRUE), simplify=TRUE, USE.NAMES= TRUE)
index_custom_Zhang <- index_custom_Zhang[index_custom_Zhang_number>2] #get rid of genesets with less than 10 genes indexed
length(index_custom_Zhang)

# perform enrichments for all contrasts
custom_enrichments_Zhang <- sapply(names(mycontrasts) , function(x) camera(y_donor, index_custom_Zhang, design_donor, mycontrasts[[x]]), simplify=FALSE, USE.NAMES= TRUE )

# add mean logfc per geneset
custom_enrichments_Zhang <- sapply(names(custom_enrichments_Zhang), function(x) 
{
  a <-  data.frame(meanlogfc = unlist(sapply(rownames(custom_enrichments_Zhang[[x]]), function(f) 
  {
    mean(tables_all_unsorted[[x]][tables_all_unsorted[[x]]$SYMBOL %in% customGeneSets_Zhang[[f]][, "SYMBOL"], "logFC"])
  }, simplify=FALSE, USE.NAMES= TRUE)))
  custom_enrichments_Zhang[[x]] <- cbind(custom_enrichments_Zhang[[x]], a )
}
,simplify=FALSE, USE.NAMES= TRUE)

custom_enrichments_Zhang <- sapply(custom_enrichments_Zhang, function(x) {x$geneset <- rownames(x); x}, simplify=FALSE, USE.NAMES= TRUE )

#show top 10 gene sets per contrast
sapply(custom_enrichments_Zhang, function(x) head(x, n = 10), simplify = FALSE)

#define boxplot function
barplot_genesets <- function(data, title){
  ggplot(data, aes(x=geneset, y=-log10(FDR), fill=meanlogfc)) +
    geom_bar(stat="identity", width =0.9, color = "black") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank() ,axis.title = element_text(size = 14), panel.background = element_rect(fill = "white", colour = "black",  size = 1/2, linetype = "solid"), axis.text = element_text(size = 16, color = "black"), legend.key.size = unit(1.1, "cm"), legend.title = element_text(size = 16),  plot.title = element_text(size= 18, hjust= 0.5 ), legend.text = element_text(size = 14)) +
    labs(x ="", y="-log10 FDR", title = title)   +
    scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))) + 
    coord_flip() 
}
#define theme
theme_marcel <- function(){ 
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size = 12), 
      legend.text = element_text(size = 12), 
      axis.title = element_text(size = 12), 
      title = element_text(size = 14)
    )
}

#define plot
make_ZhangPlots <- function(contrast = "Tex_vs_Trest")
{
  data <- custom_enrichments_Zhang[[contrast]]
  data$geneset <- factor(rownames(data), levels = rev(names(customGeneSets_Zhang)))

  data <-  data %>% dplyr::mutate(enrichment = ifelse(meanlogfc >= 0, FDR, 1)) %>% dplyr::mutate(depletion = ifelse(meanlogfc < 0, FDR,1))
  
  p <-  ggplot(data, aes(x=geneset, fill=meanlogfc)) +
    geom_bar(stat="identity", width =0.9, color = "black", aes(y= -log10(enrichment)) ) +
    geom_bar(stat="identity", width =0.9, color = "black", aes(y= log10(depletion)) ) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank() ,
          axis.title = element_text(size = 14), panel.background = element_rect(fill = "white", colour = "black",  size = 1/2, linetype = "solid"),
          axis.text = element_text(size = 16, color = "black"), legend.key.size = unit(1.1, "cm"), 
          legend.title = element_text(size = 16),  plot.title = element_text(size= 18, hjust= 0.5 ), 
          legend.text = element_text(size = 14)) +labs(x ="", y="-log10 FDR", title = paste(contrast, "Enrichment within Zhang et al 2021"))  +
    scale_y_continuous(breaks = seq(-14, 22, 2),  labels = as.character(c(seq(14, 0, -2), seq(2, 22, 2)))) +  
    scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"RdBu"))(100))) + coord_flip() 
  print(p)
}

#for every contrast plot the genes 
pdf("./output/plots/scRNAseq_Enrichment_Zhang.pdf", width=9, height = 5)
for(i in 1:6){
  make_ZhangPlots(names(custom_enrichments_Zhang)[i])
}
dev.off()

#heatmaps of genes in gene sets
pdf(file = "./output/plots/custom_enrichments_ZhangGeneset_heatmaps.pdf", width = 7, height = 8)
for(i in 1:length(customGeneSets_Zhang))
{
  gs <- names(customGeneSets_Zhang)[i]
  #plot these genes in a heatmap for each sample
  genes_of_interest <- cpms_filtered[rownames(cpms_filtered) %in% customGeneSets_Zhang[[gs]][,"SYMBOL"], ]
  topAn <- HeatmapAnnotation(
    df = data.frame(group = as.numeric(y_donor$samples[,"SampleGroup"])),  
    col = list(group = c("1" = "#F8F9F9", "2"= "#BFBFC1","3"=  "#FB0207", "4"= "#800002") ),
    annotation_legend_param = list(group = list(labels=levels(y_donor$samples$SampleGroup), title="Condition")),
    annotation_name_gp = gpar(fontsize = 0)
  )
  
  ht <- Heatmap(t(scale(t(genes_of_interest), scale=TRUE, center=TRUE)), col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)), name = "z-score log(cpm)",
                column_title =  paste(gs, "\nZhang Gene Set") , column_title_gp = gpar(fontsize=14), cluster_columns  = TRUE, cluster_rows = TRUE,  
                show_column_names = FALSE,show_row_names = TRUE,show_row_dend = FALSE, 
                row_names_gp = gpar(fontsize = 4.5, justify = "left"),top_annotation = topAn,
                heatmap_legend_param = list(labels_gp = gpar(fontsize = 12 ),
                                            title_gp = gpar(fontsize = 12 )), use_raster = FALSE)
  draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  
  write.table(genes_of_interest, file = paste0("output/lists/Zhang_Heatmaps/", gs, ".txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
}
dev.off()


#### #melt the custom_enrichment lists to perform summary plots
custom_enrichments_Zhang_melted <- sapply(names(custom_enrichments_Zhang), function(x) {custom_enrichments_Zhang[[x]]$contrast = x; custom_enrichments_Zhang[[x]]}, simplify = FALSE )
custom_enrichments_Zhang_melted <- do.call(rbind.data.frame, custom_enrichments_Zhang_melted)
rownames(custom_enrichments_Zhang_melted) <- 1:nrow(custom_enrichments_Zhang_melted)

custom_enrichments_Zhang_melted$geneset <- factor(custom_enrichments_Zhang_melted$geneset, levels =  rev(names(customGeneSets_Zhang)))
custom_enrichments_Zhang_melted$contrast <- factor(custom_enrichments_Zhang_melted$contrast, levels =  names(mycontrasts))
head(custom_enrichments_Zhang_melted)

write.table(custom_enrichments_Zhang_melted, file = "output/lists/Zhang_enrichments_list.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


# Dot plot of each gene set for the selected conditions (Figures)
pdf(file = "./output/plots/custom_genesets_Zhang_dotplot_all.pdf", width = 5, height = 4.5)
ggplot(custom_enrichments_Zhang_melted, aes(x = contrast , y = geneset, size = -log10(FDR), color = meanlogfc)) +
  geom_point() +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_size(range = c(0.3,6))+
  theme(axis.text.x =  element_text (angle = 45, hjust = 1),
        axis.text = element_text(color = "black", size = 12),
        legend.text =element_text(color = "black", size = 12)  ) 

custom_enrichments_Zhang_melted %>% filter(contrast %in% c("Ttumor_vs_Trest", "Teff_vs_Trest", "Tex_vs_Trest")) %>%
ggplot( aes(x = contrast , y = geneset, size = -log10(FDR), color = meanlogfc)) +
  geom_point() +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_size(range = c(0.3,6))+
  theme(axis.text.x =  element_text (angle = 45, hjust = 1),
        axis.text = element_text(color = "black", size = 12),
        legend.text =element_text(color = "black", size = 12)  ) 

dev.off()

###########################################################################################
# overlap analysis for hit selection (based on human TIL 2018 data)
##########################################################################################
customGeneSets <- list()
for (f in list.files("./input/custom_genesets/humanTIL2018/", pattern = "*.txt")) {
  myGeneSet <- read.table(paste0("./input/custom_genesets/humanTIL2018/", f), h=F, sep="\t")
  cat("  ", length(myGeneSet[,1]), "genes in geneset ")
  
  ## If you read the genes as symbols in the files, match them to y$genes abd use the ENTREZID column
  
  ## Add gene set to list
  clean_names <- gsub(".txt", "", f)
  clean_names <- gsub(x = clean_names, pattern = "-", replacement = "_")
  customGeneSets[[clean_names]] <- myGeneSet[,1]
  cat(clean_names, "\n")
}
length(customGeneSets)


selectedHits <- read.table("./input/selected_genes.txt", sep="\t", header= FALSE)$V1

anova_all.signifGenes[anova_all.signifGenes$SYMBOL %in% selectedHits,]

genes_signif <- anova_all.signifGenes[anova_all.signifGenes$logFC.Tex >0 |anova_all.signifGenes$logFC.Teff >0 , c("ENTREZID", "SYMBOL", "logFC.Tex", "FDR")]
dim(genes_signif)
genes_signif

selectedHits[!selectedHits %in% genes_signif$SYMBOL]

shared_matrix_genes_signif <- data.frame(sapply(names(customGeneSets),function(x) genes_signif$ENTREZID %in% customGeneSets[[x]]))
shared_matrix_genes_signif <- data.frame(genes_signif,  "nr_overlaps" = rowSums(shared_matrix_genes_signif), ifelse(shared_matrix_genes_signif=="TRUE", 1, 0))
shared_matrix_genes_signif <- shared_matrix_genes_signif[order(shared_matrix_genes_signif[,"nr_overlaps"], decreasing = TRUE),]
#shared_matrix_genes_signif <- shared_matrix_genes_signif[order(shared_matrix_genes_signif[,"nr_overlaps_ThommenPD1high_up"], decreasing = TRUE),]
head(shared_matrix_genes_signif)
combined_shared_lists_select <- cbind(shared_matrix_genes_signif, "selectedHit" = shared_matrix_genes_signif$SYMBOL %in% selectedHits)
combined_shared_lists_select <- combined_shared_lists_select[order(combined_shared_lists_select$selectedHit == TRUE, decreasing = TRUE),]
write.table(combined_shared_lists_select ,"./output/lists/shared_matrix_genes_signif_selected.txt", sep ="\t", quote=FALSE, row.names = FALSE)










