# CXCL4 Post-Integration (after Soupx,Demux, Doubletfinder) Leukocyte Subclustering and Analysis
library(Seurat)
library(ggplot2)
library(dplyr)
library(ROCR)
library(parallel)
library(ggpubr)
library(writexl)
library(readxl)
library(stringr)
library(gprofiler2)
library(tidyverse)
library(ggrepel)
library(dittoSeq)
windowsFonts("Arial" = windowsFont("Arial"))

setwd("c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/ecm/Leuko_v2/")
sample="CXCL4_Sik_Leukos_"
Samples.combined <- readRDS(file = "c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/ecm/leukocytes_cxcl4_subclusters_scvi_output_clustered.h5ad.ECM.rds") #if its an rds file

# Step 1: Subset only leukocytes clusters and remove all non-MPC clusters
Idents(Samples.combined)<- "leiden0.2"
Samples.combined = subset(Samples.combined, idents = c("Mac1", "Tcells", "Mac2", "Bcells"))

#Count number of cells
After_Subset<-length(Samples.combined@meta.data$orig.ident)
write.csv(After_Subset, file=paste0(sample,"_Cellnumber.csv"), row.names = F)

# Step 2: Renormalization and Find Variable Features
Samples.combined <- NormalizeData(Samples.combined, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Samples.combined <- FindVariableFeatures(Samples.combined, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose=FALSE)
Samples.combined <- RunUMAP(Samples.combined, reduction="X_scvi", dims=1:20)

# Save RDS
saveRDS(Samples.combined, file=paste0(sample, "03_seurat_obj_after_integration.rds"))

# Step 3: Dimplot and Compositional Analysis
# Rename and reorder Macs
Idents(Samples.combined)<-"leiden0.2"
current.cluster.ids<-levels(x=Samples.combined)
new.cluster.ids<-c("Mac", "T-Cells", "B-Cells", "Spp1 Mac ")
Samples.combined$leiden0.2<- plyr::mapvalues(x=Samples.combined$leiden0.2, from = current.cluster.ids, to=new.cluster.ids)
my_levels<- c("Mac", "Spp1 Mac ","B-Cells", "T-Cells")
Samples.combined@meta.data$leiden0.2 <- factor(x=Samples.combined@meta.data$leiden0.2, levels=my_levels)
rm(current.cluster.ids, new.cluster.ids, my_levels, After_Subset)

DimPlot(Samples.combined, reduction = "umap", group.by = "leiden0.2" , pt.size = 1, cols = c("#D9AF6B", "#d96b95","#6bd9af", "#6699cc"))
ggsave(filename=paste0(sample, "01_Dimplot_BeforeIntegration.jpeg"), width=5 , height = 4)

dittoBarPlot(Samples.combined, "leiden0.2",group.by = "condition", color.panel = c("#6bd9af","#D9AF6B", "#d96b95", "#6699cc"), legend.show = F)+ NoLegend()+ theme(
  axis.title = element_blank(),
  text = element_text(size=8, family = "Arial"),
  axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=8, family = "Arial"),
  plot.title = element_text(size=9, family="Arial", face = "plain"),
  axis.line = element_line(size = 0.2), 
  axis.ticks = element_line(size = 0.2))
ggsave(filename=paste0(sample, "02_Dittobarplot.jpeg"), width=4 , height = 4)
ggsave(filename=paste0(sample, "02_Dittobarplot.svg"), width=4, height = 4)

# Step 4: Top 5 Marker genes as Dotplot
Idents(Samples.combined)<-"leiden0.2"
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3)
write_xlsx(all.markers,path = paste0(sample,"04_Marker_AllMarkersTable.xlsx"))
top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% pull(gene)
top5 <-unique(top5)
DotPlot(Samples.combined, features = top5,dot.scale = 3)+       
  coord_flip()+ theme(axis.title = element_blank(), 
                      text = element_text(size=8, family = "Arial"),
                      axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
                      axis.text.y = element_text(size=8, family = "Arial"),
                      legend.text = element_text(size=6, family = "Arial"))
ggsave(filename = paste0(sample, "04_Marker_Dotplot.jpeg"), width = 2.6,height = 3.5, dpi = "retina", bg = "white")
ggsave(filename = paste0(sample, "04_Marker_Dotplot.svg"), width = 2.6,height = 3.5,  dpi = "retina")


# Step 5: Violin Plots for Marker genes Fn1, Spp1 and C1qa
features=c("Fn1", "Spp1", "C1qa")
Leuko_colors=c("#D9AF6B", "#d96b95", "#6bd9af","#6699cc")
for (f in features) {
     VlnPlot(Samples.combined, features = f,  group.by="leiden0.2", cols = Leuko_colors, pt.size = 0) + NoLegend()+ 
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)+theme(
        axis.title = element_blank(),
        text = element_text(size=8, family = "Arial"),
        axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
        axis.text.y = element_text(size=8, family = "Arial"),
        plot.title = element_text(size=9, family="Arial", face = "plain"),
        axis.line = element_line(size = 0.2), 
        axis.ticks = element_line(size = 0.2))
      ggsave(filename = paste(sample,"Vln_",f,".jpeg"), width=1.5 , height = 1.7)
      ggsave(filename = paste(sample,"Vln_",f,".svg"), width=1.5 , height = 1.7)
}

# Step 6: ECM Regulator Scoring

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

# ECM Regulator Scoring: Add average expression of genes in gset minus the average expression of 35 random genes 
Idents(Samples.combined)<-"leiden0.2"
ctrl_genes = 35 #important
features = matrisome_mm_genesetlist["ECM_Regulators"]
Samples.combined = AddModuleScore(object = Samples.combined, features = features, name = "ECM_Regulators", ctrl = ctrl_genes)
VlnPlot(Samples.combined, features = "ECM_Regulators1", pt.size = 0, group.by = "leiden0.2", cols = Leuko_colors) +
    geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) + NoLegend()+ theme(
      axis.title = element_blank(),
      text = element_text(size=8, family = "Arial"),
      axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
      axis.text.y = element_text(size=8, family = "Arial"),
      plot.title = element_text(size=9, family="Arial", face = "plain"),
      axis.line = element_line(size = 0.2), 
      axis.ticks = element_line(size = 0.2))
ggsave(filename = paste0(sample,"_ECMRegulators_vln_leiden02.jpeg"), width=1.8 , height = 2.1)
ggsave(filename = paste0(sample,"_ECMRegulators_vln_leiden02.svg"), width=1.8 , height = 2.1)
#ttest_calculations
Naba.df<-Samples.combined[[c("leiden0.2","condition", "ECM_Regulators1")]]
sink(paste0("ECM scoring output/",sample,"_", gset, "TTest_results.txt"))
print(paste0("t-test for ", gset, " Mac vs Spp1 Mac "))
print(t.test( Naba.df[which( Naba.df$leiden0.2=="Mac"),3], Naba.df[which( Naba.df$leiden0.2=="Spp1 Mac "),3], alternative = "two.sided"))
sink()
