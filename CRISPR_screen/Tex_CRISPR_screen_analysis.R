#This is the script to analyze the CRISPR 

####################################
# Setup
####################################
library(edgeR)
library(RColorBrewer)
library(viridis)
library(reshape2)
myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
library(ComplexHeatmap)
library(stringr)
library(circlize)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(GGally)
library(gridExtra)
library(grid)
set.seed(123)
setwd("/Volumes/CIMM$/Trefny/NGS/exhaustion_model/CRISPRScreens/Exp032_CRISPRScreen2_sequencing/FinalAnalysis/deposition/")
theme_marcel <- theme(plot.title = element_text(hjust=0.5, size = 16), axis.title = element_text(size = 16), axis.text = element_text(size = 14), legend.key.size = unit(1.1, "cm"),legend.text = element_text(size = 14), legend.title = element_text(size = 16), panel.background = element_rect(fill = "white", colour = "black",  size = 1/2, linetype = "solid"))

####################################
# Import Data
####################################
data.list <- list()
genelist <- list()
names.x <- c("positive119","positive122","positive124","positive141","positive144" )
names.y <- c("negative119","negative122","negative124","negative141","negative144" )
filepath <- c("./PinAPLpy_STARS_FDR025/Analysis/01_Alignment_Results/")
list.files(filepath)
for (x in 1:5){
  
  import <- read.csv(paste0(filepath, "ReadCounts_per_sgRNA/CD107a_positive_", x, "_GuideCounts.txt"), header = FALSE,sep="\t")
  if(x == 1){ 
    genelist[["guide"]] <- import$V1
    genelist[["gene"]] <- import$V2
  }
  data.list[[names.x[x]]] <- import$V3
  import <- read.csv(paste0(filepath, "ReadCounts_per_sgRNA/Control_", x, "_GuideCounts.txt"), header = FALSE,sep="\t")
  data.list[[names.y[x]]] <- import$V3
}
data <- matrix(unlist(data.list), nrow=length(data.list[[1]]))
data <- cbind(data.frame(unlist(genelist$guide), unlist(genelist$gene), stringsAsFactors = FALSE), data)
colnames(data) <- c("guide", "gene", names(data.list))
head(data)

write.table(data, file = "./output/lists/annotated_GuideCounts_all_samples.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
#build DGEList from the count data
# data <- data[!data$gene %in% c("LAT", "ZAP70", "LAMP1"),]
y <- DGEList(counts = data[,3:ncol(data)], genes=data[,1:2])
y <- calcNormFactors(y)
y

####################################
# Cumulative Distributions
####################################
pdf(paste0("./output/plots/cumulativeDistributions_sgRNAs.pdf"), width = 10)
y_averages <- data.frame("CD107a_positive" =  rowSums(y$counts[,grepl("positive", colnames(y$counts))]), "CD107a_negative" =  rowSums(y$counts[,grepl("negative", colnames(y$counts))]) )
distribution <- sapply(y_averages, function(x) sort(x), USE.NAMES  = TRUE)
distribution_cumsum <- apply(distribution,2, function(x) cumsum(x)/max(cumsum(x)))
#distribution_cumsum <- cbind(distribution_cumsum,"position" = seq(1, nrow(distribution_cumsum), 1) )
distribution_cumsum <- melt(distribution_cumsum)
distribution_cumsum <- data.frame("condition" = distribution_cumsum$Var2, "cumulative_fraction_Reads" = distribution_cumsum$value, "Cumulative_fraction_sgRNA" = distribution_cumsum$Var1/max(distribution_cumsum$Var1))

ggplot(distribution_cumsum, aes(x=Cumulative_fraction_sgRNA, y=cumulative_fraction_Reads)) + geom_line(aes(color = condition), size = 1) + geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5) + labs(y = "Cumulative Fraction Reads", x= "Cumulative Fraction sgRNAs")+ scale_color_brewer(palette="Dark2") + scale_x_continuous(limits = c(0, 1), expand = c(0,0)) + scale_y_continuous(limits = c(0, 1), expand = c(0,0))  + ggtitle("Cumulative Distributions sgRNAs") + theme_marcel

#make histograms of counts for each sample
counts.melted <- melt(data.frame(y$counts, "guides" = y$genes$guide))
colnames(counts.melted) <- c("guides","sample", "count" )
p <- ggplot(data= counts.melted, aes(x = count, col = sample)) + geom_freqpoly() + facet_wrap(~sample,ncol=2)
print(p)
#minimum counts
apply(y$counts, 2, function(x) min(x) )
#maximum counts
apply(y$counts, 2, function(x) max(x))
dev.off()


####################################
# PCA
####################################
CRISPR_cpms <- cpm(y) #retrieve counts per million
# define the grouping factors
condition <- factor(substring(colnames(y), 1, sapply(colnames(y), nchar)-3))
donor <- factor(substring(colnames(y), sapply(colnames(y), nchar)-2, sapply(colnames(y), nchar)))
group <- paste0(  condition, "-", donor)

y$samples <- cbind(y$samples, condition)
y$samples <- cbind(y$samples, donor)
y$samples <- cbind(y$samples, group)

par(mfrow=c(1,2))
#  MDS Plots
pdf(paste0("./output/plots/MDS_PCA.pdf"), width = 7, height= 6)
mds <- plotMDS(y, col= myPalette[as.numeric(y$samples$condition)], main="MDS Plot All Samples Exp032")
mds_plot <- data.frame("dim1" = mds$x, "dim2"= mds$y, "condition"= as.character(condition), "donor" =as.character(donor), "group" = str_replace_all(as.character(group), pattern= "-", " "))
mds_plot
ggplot(mds_plot, aes(x=dim1, y=dim2, color=condition, label=group )) + geom_point(size = 3) + geom_text_repel(show.legend = FALSE, size = 6)  + scale_color_brewer(palette = "Dark2") + labs(x = "Leading logFC dim1", y= "Leading logFC dim2") + theme_marcel + ggtitle("CRISPR Screen MDS Plot")

#  PCA Plots
pca1 <- prcomp(t(CRISPR_cpms), scale = T)
summary(pca1) ## 59% variance on PCs 1:3

n= 1
m= 2
pchs <- 21
col.v <- ifelse(condition == "positive", "#F8F9F9","#D9D9D9" )

scores <- pca1$x

scores <- data.frame(scores, "condition"= as.character(condition), "donor" =as.character(donor), "group" = str_replace_all(as.character(group), pattern= "-", " "))
scores
ggplot(scores, aes(x=PC1, y=PC2, color=condition, label=group )) + geom_point(size = 2) + geom_text_repel(show.legend = FALSE)  + labs(y= paste("PC", 2, ": ", round(summary(pca1)$importance[2,2],3)*100, "% variance explained", sep="") , x= paste("PC", 1, ": ", round(summary(pca1)$importance[2,1],3)*100, "% variance explained", sep="")) + theme_marcel

dev.off()


####################################
# Fold Change Analysis - Donor Centered
####################################

foldchanges_donors <- data.frame("guide" = y$genes$guide,
                                 "gene" = y$genes$gene,
                                 "positive119_vs_negative119" = apply(CRISPR_cpms, 1, function(x) x[1]/x[2]), 
                                 "positive122_vs_negative122"= apply(CRISPR_cpms, 1, function(x) x[3]/x[4]), 
                                 "positive124_vs_negative124" = apply(CRISPR_cpms, 1, function(x) x[5]/x[6]), 
                                 "positive141_vs_negative141" =apply(CRISPR_cpms, 1, function(x) x[7]/x[8]),
                                 "positive144_vs_negative144" =apply(CRISPR_cpms, 1, function(x) x[9]/x[10]))

foldchanges_donors.log <- data.frame(foldchanges_donors[,1:2], log(foldchanges_donors[,3:ncol(foldchanges_donors)], 2))
head(foldchanges_donors.log)
write.table(foldchanges_donors, paste0("./output/lists/foldchanges_donors.txt") ,sep='\t', quote = FALSE , row.names = TRUE, col.names = NA)
Heatmap(cor(foldchanges_donors.log[,3:6]) )
foldchanges_donors.log.melted <- melt(foldchanges_donors.log)

foldchanges_donors.log.melted$donor <- gsub(".*negative", "", foldchanges_donors.log.melted$variable)
foldchanges_donors.log.melted

#calculate the mean logfc for each gene (among its guides) for all donors
foldchanges_donors.log.melted.perGene  <- foldchanges_donors.log.melted %>%  group_by(gene, donor) %>% dplyr::summarize(meanlogfc_perDonor = mean(value))

#based on median
generanking_perGene <-  foldchanges_donors.log.melted.perGene %>% group_by(gene) %>% summarize(median =  median(meanlogfc_perDonor)) %>% arrange(desc(median))  
generanking_perGene$gene<- factor(generanking_perGene$gene, levels = as.character(generanking_perGene$gene))
foldchanges_donors.log.melted.perGene$gene <- factor(foldchanges_donors.log.melted.perGene$gene, levels = c(as.character(generanking_perGene$gene)))
foldchanges_donors.log.melted.perGene<- left_join(foldchanges_donors.log.melted.perGene,generanking_perGene, by = "gene" )

#export the data for prism, sorted based on median
write.table(foldchanges_donors.log.melted.perGene, file = "./output/lists/log2FoldChanges_perDonor.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(generanking_perGene, file = "./output/lists/generanking_perGene.txt", sep = "\t", quote = FALSE, row.names = FALSE)

