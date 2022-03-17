#Post-Integration Analysis of MPC after Cluster Exclusion
library(Seurat)
library(cowplot)
library(ggplot2)
library(clustree)
library(genesorteR, quietly = TRUE)
library(writexl)
library(dplyr)
library(ggpubr)
library(dittoSeq)
library(reshape2)
library(scico)
library(reticulate)
library(readxl)
library(stringr)
library(pals)

#Working directory
# Load MPC Sample
setwd("c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/Forte_MI_data/")
Forte_MI_dataall <- readRDS(file = "c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/Forte_MI_data/Forte_MI_data_first_clustering.rds") #if its an rds file
source(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/Volcano_plot_KH.R")
sample="MI_dataAll"

#load all samples and merge
sc.sham1= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/Sham1_out/outs/filtered_feature_bc_matrix/")
sc.sham2= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/Sham2_out/outs/filtered_feature_bc_matrix/")
sc.none= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/none_out/outs/filtered_feature_bc_matrix/")
sc.MIday1= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday1_out/outs/filtered_feature_bc_matrix/")
sc.MIday3= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday3_out/outs/filtered_feature_bc_matrix/")
sc.MIday5= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday5_out/outs/filtered_feature_bc_matrix/")
sc.MIday7= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday7_out/outs/filtered_feature_bc_matrix/")
sc.MIday14= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday14_out/outs/filtered_feature_bc_matrix/")
sc.MIday28= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday28_out/outs/filtered_feature_bc_matrix/")
sc.list = list(sc.sham1,sc.sham2,sc.MIday1,sc.MIday3,sc.MIday5,sc.MIday7,sc.MIday14,sc.MIday28)
names(sc.list) = c("sc.sham1","sc.sham2","sc.MIday1","sc.MIday3","sc.MIday5","sc.MIday7","sc.MIday14","sc.MIday28")
for (sc.sample in names(sc.list)) {sc.list[[sc.sample]]=CreateSeuratObject(counts = sc.list[[sc.sample]], project = sc.sample, min.cells = 3, min.features = 200)}
sc.none <- CreateSeuratObject(counts = sc.none, project = "sc.none", min.cells = 3, min.features = 200)
samples.merged <- merge(x = sc.none, y = sc.list , project = "Forte_MI_merged")
rm(sc.sham1,sc.sham2,sc.MIday1,sc.MIday3,sc.MIday5,sc.MIday7,sc.MIday14,sc.MIday28)

# Step 1: As described by Forte et al: Filterting, processing and initial tSNE with Clustering
samples.merged[['percent.mt']] = PercentageFeatureSet(samples.merged, pattern = '^mt-')
samples.merged =  subset(samples.merged, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA<15000)
samples.merged <- NormalizeData(samples.merged, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
samples.merged <- FindVariableFeatures(samples.merged, selection.method = 'vst', nfeatures = 2000,y.cutoff=0.5, verbose = FALSE)
samples.merged <- ScaleData(samples.merged, verbose = FALSE)
samples.merged <- RunPCA(samples.merged, verbose = FALSE)
samples.merged <- RunTSNE(samples.merged, reduction = "pca", dims = 1:24)
samples.merged <- FindNeighbors(samples.merged, reduction = "pca", dims = 1:24)
samples.merged <- FindClusters(samples.merged, resolution = 0.5)

#Select Optimal dimension for later UMAP representation
pdf(file = paste0(sample, '_02_dimension tryout.pdf'),paper = "a4",height = 23,width = 33)
dim.range = seq(5,50,by=5)
dim.range.umaps = vector("list")
dim.range.umaps[[1]] <- ElbowPlot(Forte_MI_dataall,50)
for (dim in dim.range) {
  message(paste0("Testing step for dim: ",dim))
  Forte_MI_dataall <- RunUMAP(Forte_MI_dataall, dims = 1:dim)
  dim.range.umaps[[paste0("dim:",dim)]]<-DimPlot(Forte_MI_dataall, group.by = "seurat_clusters") + 
    NoLegend() + labs(title=paste0("dim: 1:",dim))+theme(axis.title.y = element_blank(),axis.text = element_text(size = 8),axis.title = element_text(size = 10))
}
print(ggarrange(plotlist = dim.range.umaps, widths = 20, heights = 30,ncol = 3,nrow = 3))
dev.off()

# Step 2: Annotation of clusters using top 10 Marker Genes
all.markers = FindAllMarkers(Forte_MI_dataall, test.use = "MAST",min.pct = 0.3,assay = "RNA")
write_xlsx(all.markers, path=paste0(sample, "_10_allmarkers.xlsx"), col_names = TRUE)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene)
top10 <-unique(top10)
DotPlot(Forte_MI_dataall, features = top10,dot.scale = 10, group.by = "clusternames") + coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(filename = paste0(sample,"_11_top10_featureDOTplot.jpeg"), width=10 , height = 23)
# Annotate Clusters, remove Doublets
current.cluster.ids = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25")
new.cluster.ids = c("Res Mac", "EC #1", "Fibro #1", "Myofibro", "Matrifibro", "Fibro #2", "Ly6c2hi Mono", "Granulocytes", "Progen. Fibro", "Arg1 Mac", "B Cells", "EndD Fibro", "Fibro #3", "Prolif. Fibro", "DC", "NK_T Cells", "SMC_Peri", "EC #2", "EC #3", "Prolif. Mac","IFN Fibro", "Lymphatic EC", "Fibro/Mac/EC Doublets", "Schwann-Cells", "B Cell/Fibro Doublets", "B-Cell/EC Doublets")
Forte_MI_dataall$Annotation_Level_3<- plyr::mapvalues(x=Forte_MI_dataall$seurat_clusters, from = current.cluster.ids, to=new.cluster.ids)
Forte_MI_dataall = subset(Forte_MI_dataall, idents = c(0,1,2,3,4,5,6,7,9,8,10,11,12,13,14,15,16,17,18,19,20,21,23))

# Step 3: Dimplot of annotated clusters
Forte_MI_dataall <- RunUMAP(Forte_MI_dataall, dims=1:30)
DimPlot(Forte_MI_dataall, reduction = "umap", group.by = "Annotation_Level_3" , pt.size = 0.8, label = TRUE, repel = TRUE,label.size = 5)
ggsave(filename=paste0(sample, "_13_ClusterNames_Level3.jpeg"),width=9, height= 8)

## Feature and Count plots
FeaturePlot(Forte_MI_dataall, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_08_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Forte_MI_dataall,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_09_nfeature_ncount_vln.jpeg"), width=10 , height =5)

#Count number of cells
After_Subset<-length(Forte_MI_dataall@meta.data$orig.ident)
write.csv(After_Subset, file=paste0(sample,"_Cellnumber.csv"), row.names = F)

# Step 4: Top3 Marker Genes as Dotplot
Idents(Forte_MI_dataall)<-"Annotation_Level_3"
all.markers = FindAllMarkers(Forte_MI_dataall, test.use = "MAST",min.pct = 0.3,assay = "RNA")
write_xlsx(all.markers, path=paste0(sample, "_15_allmarkers.xlsx"), col_names = TRUE)
top3 <- all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top3 <-unique(top3)
DotPlot(Forte_MI_dataall, features = top3,dot.scale = 2, group.by = "Annotation_Level_3") + coord_flip() + theme(
  axis.title = element_blank(),
  text = element_text(size=6, family = "Arial"),
  axis.text.x = element_text(size=6, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=6, family = "Arial"))
ggsave(filename = paste0(sample,"15_top3_featureDotplot.jpeg"), width=4 , height = 6, bg = "white")
ggsave(filename = paste0(sample,"15_top3_featureDotplot.svg"), width= 4, height = 6)

# Step 5: Save Seurat Object as well as subsetted immune cells for downstream analysis
saveRDS(Forte_MI_dataall, file = paste0(sample, "_13_SeuratObjt_AfterFirstAnalysis.RDS"))
Idents(Forte_MI_dataall) <- "seurat_clusters"
Forte_MI_MPC <- subset(Forte_MI_dataall, idents=c("0","6","9", "14"))
saveRDS(Forte_MI_MPC, file=paste0(sample, "_14_SeuratObj_ImmuneCells.RDS"))
