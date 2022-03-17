# CXCL4 Post-Integration Corrrelation analysis of annotated leukocytes with Forte MI Leukocytes
windowsFonts("Arial" = windowsFont("Arial"))
library(Seurat)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(readxl)
options(future.globals.maxSize= 1000*1024^2)

setwd("c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/ecm/Leuko_v2/")
sample="CXCL4_SikVsForte_Cor"
sc_snRNA <- readRDS(file = "c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/ecm/Leuko_v2/CXCL4_Sik_Leukos_03_seurat_obj_after_integration.rds") #if its an rds file
sc_Forte <- readRDS("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/Forte_MI_data/ImmuneCellAnalysis/MI_ImmuneCell_13_SeuratObjt_AfterFirstAnalysis.RDS")
sc_Forte@assays$ECMRegulators=NULL # Remove calculated ECM Regulator scores as they are not needed)

# Subset Day 28 of Forte MI Dataset as we only have day 28 for snRNA data to enable better comparison
Idents(sc_Forte) <- "timepoint"
sc_Forte <- subset(sc_Forte, idents="MI_Day28")

# Remove assays except for RNA
sc_Forte@assays$integrated=NULL # Remove calculated ECM Regulator scores as they are not needed)

# Rename assay for merging and place idents under Annotation_Level_1
sc_snRNA<-RenameAssays(object = sc_snRNA, "originalexp" = 'RNA')
sc_snRNA$Annotation_Level_1 <- sc_snRNA$leiden0.2

# Merge datasets processing
sc_merged <- merge(x = sc_snRNA, y = sc_Forte, project = sample)
sc_merged <- sc_merged %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose=FALSE) %>%
  RunPCA(pc.genes = sc_merged@var.genes, npcs = 20, verbose = FALSE)
  
# Pull top 3 Marker genes for annotated leukocytes from Forte et al and annotated leukocytes from snRNA dataset
Idents(sc_Forte) <- "Annotation_Level_1"
all.markers.forte<-FindAllMarkers(sc_Forte, test.use = "MAST",min.pct = 0.3, assay = "RNA")
top3.forte <- all.markers.forte %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)

all.markers.snRNA = read_excel("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/ecm/Leuko_v2/CXCL4_Sik_Leukos_04_Marker_AllMarkersTable.xlsx")
top3.snRNA <- all.markers.snRNA %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top6<- c(top3.forte, top3.snRNA)
top6 <-unique(top6) # perform correlation analysis on a total of 35 differentially expressed genes across clusters

#-----pearson correlation of subclusters based on the marker genes identified using FindAllMarkers functions for both datasets----------------------
Idents(sc_merged)="Annotation_Level_1"
avgexp = AverageExpression(sc_merged, return.seurat = T)
expression.values <-FetchData(avgexp,vars =top6, slot = "scale.data")
expression.values <-t(expression.values)
avgexp.cor = cor(x = expression.values[,1:4],y = expression.values[,5:12])
avgexp.cor.heatmap <- pheatmap(avgexp.cor,main = "pearson_correlation", cluster_rows = T, cluster_cols = T, angle_col = 45,border_color = 0,treeheight_row = 0, treeheight_col = 0, fontsize = 14, fontsize_row = 14)
jpeg(filename = paste0(sample,"_pearson_correlation_of_subclusters.jpeg"), width=1200 , height =1000,quality = 100,res = 300)
print(avgexp.cor.heatmap)
dev.off()
setEPS()
postscript(paste0(sample,"_pearson_correlation_of_subclusters.eps"))
print(avgexp.cor.heatmap)
dev.off()
