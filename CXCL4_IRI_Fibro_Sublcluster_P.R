# CXCL4 Post-Integration (after Soupx,Demux, Doubletfinder) Fibroblast Subclustering and Analysis
library(Seurat)
library(ggplot2)
library(zellkonverter)
library(dplyr)
library(ggpubr)
library(writexl)
library(scales)
library(gprofiler2)
library(tidyverse)
library(ggrepel)
require(svglite)
library(dittoSeq)
library(readxl)
library(ggrastr)
windowsFonts("Arial" = windowsFont("Arial"))

setwd("c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/ecm/Fibro/")
sample="CXCL4_scvi_Fibro_"

# Convert of subsetted fibro.h5ad to seurat
sce.h5ad <- readH5AD("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/file_from_sikander/Fibroblasts_cxcl4_subclusters_scvi_output_clustered.h5ad")
Samples.combined <- as.Seurat(sce.h5ad, counts = "X", data = "X")
rm(sce.h5ad)

# Step 1: Renormalization and Find Variable Features
Samples.combined <- NormalizeData(Samples.combined, normalization.method = "LogNormalize", scale.factor = 10000)
Samples.combined <- FindVariableFeatures(Samples.combined, selection.method = "vst", nfeatures = 2000)
Samples.combined <- ScaleData(Samples.combined, verbose=FALSE)
DimPlot(Samples.combined, reduction = "X_umap", group.by = "leiden0.3", label = TRUE)
ggsave(filename = paste0(sample, "00_leiden03_BeforeSubset.jpeg"), width=6, height=5)

# Step 2: Annotate leiden0.3 clusters using the top10 marker genes to subset fibroblasts and remove doublet contamintating clusters
Idents(Samples.combined)<-"leiden0.3"
# Top 10 marker genes
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3)
write_xlsx(all.markers,path = paste0(sample,"00_leiden03_Marker_AllMarkersTable_BeforeSubset.xlsx"))  
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene)
top10 <-unique(top10)
DotPlot(Samples.combined, features = top10,dot.scale = 5)+
  coord_flip()+ theme(axis.title = element_blank(),axis.text.y=element_text(size = 6))
ggsave(filename = paste0(sample, "00_leiden03_Marker_Dotplot_BeforeSubset.svg"), width = 5,height = 10, dpi = "retina")
# Annotation
Idents(Samples.combined) <- "leiden0.3"
new.cluster.names <- c("Fibro 1", "DCT","PT", "Endo", "TAL", "Fibro 2", "mtGene Cluster", "Fibro 3", "DTL", "IC", "Fibro 4", "Podo")
old.cluster.names <- levels(Samples.combined)
Samples.combined$Annotation<- plyr::mapvalues(x=Samples.combined$leiden0.3, from = old.cluster.names  , to=new.cluster.names)
Idents(Samples.combined) <- "Annotation"
Samples.combined <- subset(Samples.combined, idents=c("Fibro 1", "Fibro 2", "Fibro 3", "Fibro 4"))
Samples.combined$Annotation<-factor(x=Samples.combined@meta.data$Annotation)
rm(new.cluster.names, old.cluster.names, top10, all.markers)

#Count number of cells
After_Subset<-length(Samples.combined@meta.data$orig.ident)
write.csv(After_Subset, file=paste0(sample,"_Cellnumber.csv"), row.names = F)
rm(After_Subset)

# Adjust condition names
Idents(Samples.combined) <- "condition"
new.condition.names <- c("Cxcl4-/- Sham", "Cxcl4-/- IRI", "WT Sham", "WT IRI")
old.contidtion.names <- levels(Samples.combined)
Samples.combined$condition<- plyr::mapvalues(x=Samples.combined$condition, from = old.contidtion.names, to=new.condition.names)
Samples.combined$condition<-factor(x=Samples.combined@meta.data$condition, levels=c("WT Sham", "Cxcl4-/- Sham", "WT IRI", "Cxcl4-/- IRI"))
rm(new.condition.names, old.contidtion.names)

# Step 3: Dimplot and Dittobarplot for Composition Analysis
# Dimplot
DimPlot(Samples.combined, reduction = "X_umap", group.by = "Annotation" , pt.size = 1, cols = c("#1D6996","#0F8554","#E17C05","#6F4070"))
ggsave(filename=paste0(sample, "01_Dimplot.jpeg"), width=5 , height = 4)

# Composition Bar Plot
# create color palette from UMAP
identities <- levels(Samples.combined@meta.data$Annotation)
my_color_palette <- hue_pal()(length(identities))
# Dittobarplot for Compositional Analysis
dittoBarPlot(Samples.combined, "Annotation", group.by = "condition", color.panel = c("#1D6996","#0F8554","#E17C05","#6F4070"), x.reorder = c(4,2,3,1), legend.show = F)+ NoLegend()+ theme(
  axis.title = element_blank(),
  text = element_text(size=8, family = "Arial"),
  axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=8, family = "Arial"),
  plot.title = element_text(size=9, family="Arial", face = "plain"),
  axis.line = element_line(size = 0.2), 
  axis.ticks = element_line(size = 0.2))
ggsave(filename=paste0(sample, "02_Dittobarplot.jpeg"), width=2 , height = 2)
ggsave(filename=paste0(sample, "02_Dittobarplot.svg"), width=2 , height = 2)

# Save RDS
saveRDS(Samples.combined, file=paste0(sample, "03_seurat_obj_after_clusterexclusion.rds"))

# Step 4: Dotplot of top3 marker genes Marker gene calculations
Idents(Samples.combined)<-"Annotation"
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3)
top3 <- all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top3 <-unique(top3)
DotPlot(Samples.combined, features = top3,dot.scale = 4)  + scale_y_discrete(limits=c("Fibro 4","Fibro 3", "Fibro 2", "Fibro 1")) + theme(
          axis.title = element_blank(),
          text = element_text(size=8, family = "Arial"),
          axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
          axis.text.y = element_text(size=8, family = "Arial"),
          legend.text = element_text(size=6, family = "Arial"))
ggsave(filename = paste0(sample, "04_Marker_Dotplot_t.jpeg"), width = 3.8,height = 1.7, dpi = "retina", bg = "white")
ggsave(filename = paste0(sample, "04_Marker_Dotplot_t.svg"),  width = 3.8,height = 1.7, dpi = "retina")
rm(top3, all.markers)

# Step 5: Core Matrisome scoring of Fibroblast Subsets
# Create new metadata column which combines cluster + Surgery for later ECM Analysis of KO v WT
Idents(Samples.combined)<-"condition"
current.cluster.ids<-levels(x=Samples.combined)
new.cluster.ids<-c("Sham", "Sham", "IRI", "IRI")
Samples.combined$surgery<- plyr::mapvalues(x=Samples.combined$condition, from = current.cluster.ids, to=new.cluster.ids)
new.cluster.ids<-c("WT", "Cxcl4-/-", "WT", "Cxcl4-/-")
Samples.combined$genotype<- plyr::mapvalues(x=Samples.combined$condition, from = current.cluster.ids, to=new.cluster.ids)
Samples.combined$Annotation_surgery <- paste(Samples.combined$Annotation, Samples.combined$surgery, sep = " ")
my_levels <- c("Fibro 1 Sham", "Fibro 1 IRI","Fibro 2 Sham", "Fibro 2 IRI", "Fibro 3 Sham", "Fibro 3 IRI", "Fibro 4 Sham", "Fibro 4 IRI")
Samples.combined@meta.data$Annotation_surgery <- factor(x=Samples.combined@meta.data$Annotation_surgery, levels=my_levels)
my_levels<- c("WT", "Cxcl4-/-")
Samples.combined@meta.data$genotype <- factor(x=Samples.combined@meta.data$genotype, levels=my_levels)

# Process mouse matrisome genes from #http://matrisomeproject.mit.edu/
matrisome_mm_masterlist <- as.data.frame(read_excel("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/matrisome_mm_masterlist.xls"))
matrisome_mm_masterlist<-matrisome_mm_masterlist[c("Division","Category","Gene Symbol")]
matrisome_mm_masterlist$Division=gsub(pattern = " ",replacement = "_",x = matrisome_mm_masterlist$Division)
matrisome_mm_masterlist$Division=gsub(pattern = "-",replacement = "_",x = matrisome_mm_masterlist$Division)
matrisome_mm_masterlist$Category=gsub(pattern = " ",replacement = "_",x = matrisome_mm_masterlist$Category)
matrisome_mm_masterlist$Category=gsub(pattern = "-",replacement = "_",x = matrisome_mm_masterlist$Category)
matrisome_mm_masterlist.2<-matrisome_mm_masterlist
matrisome_mm_masterlist.2$Category[matrisome_mm_masterlist$Division=="Matrisome_associated"] = "Matrisome_associated"
matrisome_mm_masterlist.2$Category[matrisome_mm_masterlist$Division=="Core_matrisome"] = "Core_matrisome"
matrisome_mm_masterlist<-rbind(matrisome_mm_masterlist,matrisome_mm_masterlist.2)
matrisome_mm_masterlist<-matrisome_mm_masterlist[matrisome_mm_masterlist$Division!="Retired",]
rm(matrisome_mm_masterlist.2)
matrisome_mm_genesetlist = list()
for (geneset in unique(matrisome_mm_masterlist$Category)) {
  matrisome_mm_genesetlist[[geneset]] = matrisome_mm_masterlist$`Gene Symbol`[matrisome_mm_masterlist$Category==geneset]
}

# ECM scoring for the seurat object: Add average expression of genes in gset minus the average expression of 35 random genes 
ctrl_genes = 35 #important
Idents(Samples.combined) <- "Annotation_surgery"
gset="Core_matrisome"
features = matrisome_mm_genesetlist[gset]
gset.name= paste0(gset,"1")
Samples.combined = AddModuleScore(object = Samples.combined, features = features, name = gset, ctrl = ctrl_genes)
VlnPlot(Samples.combined, features = paste0(gset, '1'), pt.size = 0,group.by = "Annotation_surgery", split.by = "genotype",split.plot= T,cols = c('#00C5C0','#FA756C')) +grids(size = 0.2)+
    geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) + theme(
      axis.title.x = element_blank(),
      axis.title   = element_blank(),
      text         = element_text(size=8, family = "Arial"),
      axis.text.x  = element_text(size=8, family = "Arial", angle=45, hjust=1),
      axis.text.y  = element_text(size=8, family = "Arial"),
      plot.title   = element_text(size=9, family="Arial", face = "plain"),
      axis.line    = element_line(size = 0.2), 
      axis.ticks   = element_line(size = 0.2),
      legend.title = element_text(size=8, family = "Arial"),
      legend.text  = element_text(size=8, family = "Arial"))
ggsave(filename = paste0("ECM scoring output/",sample,"_",gset,"_vln.jpeg"), width=4, height = 2.2)
ggsave(filename = paste0("ECM scoring output/svg/",sample,"_",gset,"_vln.svg"), width=3.5 , height = 2.2)
  
#ttest calculations
Naba.df<-Samples.combined[[c("Annotation","condition", gset.name)]]
sink(paste0("ECM scoring output/",sample,"_", gset, "TTest_results.txt"))
print(paste0("t-test for ", gset, " Fibroblast 1 IRI WT vs KO"))
print(t.test(Naba.df[which((Naba.df$Annotation=="Fibro 1")&(Naba.df$condition=="WT IRI")),3],Naba.df[which((Naba.df$Annotation=="Fibro 1")&(Naba.df$condition=="Cxcl4-/- IRI")),3], alternative = "two.sided"))
print(paste0("t-test for ", gset, " Fibroblast 2 IRI WT vs KO"))
print(t.test(Naba.df[which((Naba.df$Annotation=="Fibro 2")&(Naba.df$condition=="WT IRI")),3],Naba.df[which((Naba.df$Annotation=="Fibro 2")&(Naba.df$condition=="Cxcl4-/- IRI")),3], alternative = "two.sided"))
print(paste0("t-test for ", gset, " Fibroblast 3 IRI WT vs KO"))
print(t.test(Naba.df[which((Naba.df$Annotation=="Fibro 3")&(Naba.df$condition=="WT IRI")),3],Naba.df[which((Naba.df$Annotation=="Fibro 3")&(Naba.df$condition=="Cxcl4-/- IRI")),3], alternative = "two.sided"))
print(paste0("t-test for ", gset, " Fibroblast 4 IRI WT vs KO"))
print(t.test(Naba.df[which((Naba.df$Annotation=="Fibro 4")&(Naba.df$condition=="WT IRI")),3],Naba.df[which((Naba.df$Annotation=="Fibro 4")&(Naba.df$condition=="Cxcl4-/- IRI")),3], alternative = "two.sided"))
print(paste0("t-test for ", gset, " Fibroblast (all) IRI WT vs KO"))
print(t.test(Naba.df[which(Naba.df$condition=="WT IRI"),3],Naba.df[which(Naba.df$condition=="Cxcl4-/- IRI"),3], alternative = "two.sided"))
sink()

