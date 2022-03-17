# Adjust within th Script: 1 - working directory, 2 - sample name, 3 - Directory of read 10x data function, 4 - resolution range, 5- timepoint & mouse#, 6 - Quality filter (cutoff feature_RNA, percent.mt), 6 -  dimensions (based on Elbow Plot)
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(clustree)
library(genesorteR, quietly = TRUE)
library(cowplot)
library(ggpubr)
library(writexl)
library(clustree)
library(dittoSeq)
library(harmony)
library(scales)
library(readxl)
windowsFonts("Arial" = windowsFont("Arial"))

#Step 1: Set working directory, and name of the dataset
setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/AllIntegrated/MPC/")
sample = 'Rao_MPC'

sc <- readRDS(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/AllIntegrated/Rao_08_seurat_obj_MPC.rds")

# Normalize, find variable features & scale
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sc <- FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
sc <- ScaleData(sc, verbose=FALSE)

############################## -- HARMONY INTEGRATION --#######################################
# Dimplot before Integration
sc <- RunPCA(sc, ncps=30, verbose=FALSE)
sc <- RunUMAP(sc, reduction="pca", dims=1:30)
DimPlot(sc, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, ncol=4) + NoLegend()
ggsave(filename = paste0(sample, "_01_Dimplot_beforeReInteg.jpeg"), width= 8, height = 8)

# Harmony-Integration
sc <- sc %>% RunHarmony("orig.ident", plot_convergence=TRUE)
ggsave(filename = paste0(sample, "_01_Harmony.jpeg"), width= 4, height =4)
ggsave(filename = paste0(sample, "_01_Harmony.pdf"), width= 4, height =4)

# Dimplot After Integration
sc <- RunUMAP(sc, reduction="harmony", dims=1:30)
sc <- FindNeighbors(sc, reduction = 'harmony', dims = 1:30, verbose = FALSE)
DimPlot(sc, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, ncol=2) + NoLegend()
ggsave(filename = paste0(sample, "_01_Dimplot_afterReInteg.jpeg"), width= 8, height = 8)
 
#save 4.2: Save integrated dataset before clustering
saveRDS(sc, file = paste0(sample, "_02_seurat_obj_after_Int.rds"))  

# Step 6: clusters and feature DotPlot for multiple resolutions 
FeaturePlot(sc, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_03_featureplot.jpeg"), width=20 , height = 10)

VlnPlot(sc,features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = "orig.ident",pt.size = 0) & NoLegend() & theme(
  axis.title = element_blank(),
  text = element_text(size=8, family = "Arial"),
  axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=6, family = "Arial"),
  plot.title = element_text(size=9, family="Arial", face = "plain"),
  axis.line = element_line(size = 0.2), 
  axis.ticks = element_line(size = 0.2))
ggsave(filename = paste0(sample,"_03_Gene&RNA count Filtered.jpeg"), width=8 , height = 2)
ggsave(filename = paste0(sample,"_03_Gene&RNA count Filtered.svg"), width=8 , height = 2)

# trying different resolutions
message(paste0("start tryout for resolutions sample: ",sample))
pdf(file = paste0(sample, '_06_dimension tryout.pdf'),paper = "a4",height = 23,width = 33)
# PCA testing
message("start testing for different dimensions of the umap. Deflaut for later will be 1:30")
dim.range = seq(5,30, by=5)
dim.range.umaps = vector("list")
dim.range.umaps[[1]] <- ElbowPlot(sc,50)
for (dim in dim.range) {
  message(paste0("Testing step for dim: ",dim))
  sc <- RunUMAP(sc, reduction = "harmony", dims = 1:dim)
  dim.range.umaps[[paste0("dim:",dim)]]<-DimPlot(sc) + 
    NoLegend() + labs(title=paste0("dim: 1:",dim))+theme(axis.title.y = element_blank(),axis.text = element_text(size = 8),axis.title = element_text(size = 10))
}
print(ggarrange(plotlist = dim.range.umaps, widths = 20, heights = 30,ncol = 3,nrow = 3))
dev.off()
message(paste0("DONE with dimension tryout of ",sample))

# choose dimension of 1:20
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:20)
sc <- FindNeighbors(sc, reduction = 'harmony', dims = 1:20, verbose = FALSE)

#Clustertree
sc.int=sc
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.5, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_06_cluster_resolution_tree.jpeg"), width=10 , height = 10)

#clusters and featureDotPlot for multiple resolution
res.range = seq(0.1,0.5,by=0.2)
message(paste0("start loop for res sample: ",sample))
sc.res <- sc
for (res in res.range) {
  pdf(file = paste0(sample, "_06_resolution_", res, "_tryout.pdf"),paper = "a4",height = 23,width = 33)
  sc.res <- FindClusters(sc.res, resolution = res)
  #plot with res on umap
  p1 <- DimPlot(sc.res, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) + 
    NoLegend() + labs(title=paste0("res: ",res))+theme(axis.title.y = element_blank(),axis.text = element_text(size = 8),axis.title = element_text(size = 10))
  #plot nfeatures and ncount per cluster to identify low quality clusters
  v2 <- VlnPlot(sc.res,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0.1) +NoLegend() +theme(axis.text.x = element_text(angle = 0,hjust = 0),axis.text = element_text(size = 8),axis.title = element_text(size = 8),plot.title = element_text(size = 10))
  v3 <- VlnPlot(sc.res,group.by = "seurat_clusters", features =c('nCount_RNA'), pt.size = 0.1) +NoLegend() +theme(axis.text.x = element_text(angle = 0,hjust = 0),axis.text = element_text(size = 8),axis.title = element_text(size = 8),plot.title = element_text(size = 10))
  #top10 genes per cluster for res
  all.markers = FindAllMarkers(sc.res, only.pos = TRUE,min.pct = 0.25,assay = "RNA")
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene)
  top10 <-unique(top10)
  dp1 <- DotPlot(sc.res, features = top10,dot.scale = 5,assay = "RNA")+       
    coord_flip() + NoLegend()+ theme(axis.title = element_blank(),axis.text.y=element_text(size = 7))
  print(ggarrange(plotlist = list(ggarrange(p1,v2,v3,ncol=1,nrow = 3),dp1),ncol = 2,widths = c(1.5,2)))
  dev.off()
}
message(paste0("DONE with sample: ",sample))
rm(sc.res)

# Select resolution 0.3 and remove NK Cells
sc <- FindClusters(sc, resolution=0.3)
sc <- subset(sc, idents = 6, invert=TRUE)

# Re-run UMAP, FindNeighbors and Dimplot
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:20)
sc <- FindNeighbors(sc, reduction = 'harmony', dims = 1:20, verbose = FALSE)
identities <- levels(sc@active.ident)
my_color_palette <- c("#00C094","#ff61c3","#F8766D","#00B6EB", "#a58aff","#C49A00") # handpicked to fit/differ to huCKD MPC color_palette
DimPlot(sc, group.by = "seurat_clusters", cols = my_color_palette)
ggsave(filename=paste0(sample, "_07_Dimplot_AfterIntegration.jpeg"), width=4.5 , height = 4)

saveRDS(sc, file = paste0(sample, "_07_seurat_obj_after_Clusterexclusion.rds")) 

# Compositional Analysis heart-failure with log2FC plots
Idents(sc)<- "seurat_clusters"
identities <- levels(sc@active.ident)
my_color_palette <- hue_pal()(length(identities))
current.cluster.ids<-c("Healthy", "DCM", "ICM")
new.cluster.ids<-c("Healthy", "HeartFailure", "HeartFailure")
sc$heartfailure<- plyr::mapvalues(x=sc$disease, from = current.cluster.ids, to=new.cluster.ids)
dittoBarPlot(sc, var = "seurat_clusters", group.by = "heartfailure", scale = "percent")
dittobar_numbers=dittoBarPlot(sc, var="seurat_clusters", group.by="heartfailure", data.out = TRUE)
Cell_L1<-as.data.frame(dittobar_numbers$data)
Cell_L1_Healthy<- filter(Cell_L1, grouping=="Healthy") %>% as.data.frame() %>% select(label, percent) %>% rename(Healthy="percent")
Cell_L1_HF<-  filter(Cell_L1, grouping=="HeartFailure") %>% as.data.frame() %>% select(label, percent)  %>% rename(HeartFailure="percent")
Cell_L1 = merge(Cell_L1_Healthy,Cell_L1_HF, by = "label")
Cell_L1$HF_Log2FC <- log2(Cell_L1$HeartFailure/Cell_L1$Healthy)
Cell_L1$label <- factor(Cell_L1$label, levels = Cell_L1$label[order(Cell_L1$HF_Log2FC)])
ggplot(Cell_L1, aes(x=label, y=HF_Log2FC)) + geom_bar(stat = "identity", fill=my_color_palette) + coord_flip()
ggsave(paste0(sample, '_08_dittobar_HF.jpeg'), width =1.5, height = 3)
ggsave(paste0(sample, '_08_dittobar_HF.svg'), width =1.5, height = 3)

# Dotplot of top 3 Marker genes
all.markers = FindAllMarkers(sc, test.use = "MAST",min.pct = 0.3)
write_xlsx(all.markers,path = paste0(sample,"_09_Marker_AllMarkersTable.xlsx"))
top3 <- all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top3 <-unique(top3)
DotPlot(sc, features = top3,dot.scale = 3) + theme(axis.title = element_blank(), 
                      text = element_text(size=8, family = "Arial"),
                      axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
                      axis.text.y = element_text(size=8, family = "Arial"),
                      legend.text = element_text(size=6, family = "Arial"))
ggsave(filename = paste0(sample, "_09_Marker_Dotplot.jpeg"), width = 4.5, height = 1.8, dpi = "retina", bg = "white")
ggsave(filename = paste0(sample, "_09_Marker_Dotplot.svg"), width = 4.5,height = 1.8,  dpi = "retina")

# Calculate ECM Regulator Score 
# Process mouse matrisome genes from #http://matrisomeproject.mit.edu/
matrisome_mm_masterlist <- as.data.frame(read_excel("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/matrisome_hs_masterlist.xls"))
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
# ECM Regulator scoring for the seurat object
DefaultAssay(sc)="RNA"
ctrl_genes = 50 #important
gset="ECM_Regulators"
features = matrisome_mm_genesetlist[gset]
sc = AddModuleScore(object = sc, features = features, name = gset, ctrl = ctrl_genes)
FeaturePlot(sc, features = paste0(gset, '1'),pt.size = 1.5,order = T,label = FALSE) + 
          scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', midpoint = 0, limits = c(-1,1))
  ggsave(filename = paste0(sample,"_",gset,"_f.jpeg"), width=10 , height = 10)
  ggsave(filename = paste0(sample,"_",gset,"_f.svg"), width=10 , height = 10)
VlnPlot(sc, features = paste0(gset, '1'),  group.by="seurat_clusters",pt.size = 0, sort = "decreasing")+NoLegend() + coord_flip() +geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)& theme(
  axis.title = element_blank(),
  text = element_text(size=8, family = "Arial"),
  axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=6, family = "Arial"),
  plot.title = element_text(size=9, family="Arial", face = "plain"),
  axis.line = element_line(size = 0.2), 
  axis.ticks = element_line(size = 0.2))
ggsave(filename = paste(sample,"_ECMRegulators_V.jpeg"), width=1.6, height = 2.3)
ggsave(filename = paste(sample,"_ECMRegulators_V.svg"), width=1.6, height = 2.3)

#ttest_calculations for ECM Regulator Score & Pf4 for SPP1 Mac (Cluster 4) vs Res Mac (Cluster 0)
Naba.df<-FetchData(object = sc, vars=c("seurat_clusters", "ECM_Regulators1"))
sink(paste0("",sample,"_11_TTest_Results.txt"))
print(paste0("t-test for  ECM Regulator Spp1 Mac WT IRI vs Res Mac"))
print(t.test(Naba.df[which(Naba.df$seurat_clusters =="4"),2], Naba.df[which(Naba.df$seurat_clusters=="0"),2], alternative = "two.sided"))
sink()

# Featureplots
features_all=c("FN1", "SPP1", "APOE", "TREM2", "CD9")
dir.create("FeatureandVlnPlots")
for (f in features_all) {
      FeaturePlot(sc, features = f, reduction="umap",min.cutoff = "q9", order = TRUE,pt.size = 0.7) 
      ggsave(filename = paste("FeatureandVlnPlots/", sample,"_",f,"_f.jpeg"), width=4.5 , height = 4)
      VlnPlot(sc, features = f,  group.by="seurat_clusters", pt.size = 0) & NoLegend() & theme(
        axis.title = element_blank(),
        text = element_text(size=8, family = "Arial"),
        axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
        axis.text.y = element_text(size=6, family = "Arial"),
        plot.title = element_text(size=9, family="Arial", face = "plain"),
        axis.line = element_line(size = 0.2), 
        axis.ticks = element_line(size = 0.2))
      ggsave(filename = paste("FeatureandVlnPlots/",sample,"_",f,"_Vln.jpeg"), width=2.5 , height = 1.7)
      ggsave(filename = paste("FeatureandVlnPlots/",sample,"_",f,"_Vln.svg"), width=2.5 , height = 1.7)}
