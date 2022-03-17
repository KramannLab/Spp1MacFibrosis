# Analysis of huCD10 negative Dataset from Kuppet et al., nature for SPP1 MPC
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
library(harmony)
library(scales)
library(tibble)
library(readxl)

#Working directory and sample name
setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/huCKD_Nature_Decode_Myofibroblast_Dataset/Cxcl4")
sample="huCKD_CD10Neg_MPC_"

#Read in CD10Negative_huCKD Dataset
Samples.combined <- readRDS("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/huCKD_Nature_Decode_Myofibroblast_Dataset/Human_CD10negative.rds")
DimPlot(Samples.combined, group.by="Annotation.Level.2",pt.size = 0.8)
ggsave(filename = "huCKD_Cd10Neg_AllCells_01_Dimplot.jpeg", width=10 , height = 5)
DimPlot(Samples.combined, group.by="Annotation.Level.2",pt.size = 0.8, label = TRUE)
ggsave(filename = "huCKD_Cd10Neg_AllCells_01_Dimplot_l.jpeg", width=20 , height = 10)

#Step 1: Subset mononuclear phagocytes (MPC)
Idents(Samples.combined)<-"Annotation.Level.2"
Samples.combined <- subset(Samples.combined,idents = c("Dendritic Cells", "Monocytes", "Macrophages"))

#Step 2: Re-Integration using Harmony
Samples.combined <- NormalizeData(Samples.combined, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Samples.combined <- FindVariableFeatures(Samples.combined, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose=FALSE)
Samples.combined <- RunPCA(Samples.combined, ncps=30, verbose=FALSE)
Samples.combined <- RunUMAP(Samples.combined, reduction="pca", dims=1:30)
DimPlot(Samples.combined, reduction = "umap", group.by="Patient.ID",pt.size = 0.8, ncol=2) # Dimplot before Integration
Samples.combined <- Samples.combined %>%  RunHarmony("Patient.ID", plot_convergence=TRUE) # Run Harmony

#Step 3: Find optimal UMAP Dimension
pdf(file = paste0(sample, '_01_dimension tryout.pdf'),paper = "a4",height = 23,width = 33)
dim.range = seq(5,50,by=5)
dim.range.umaps = vector("list")
dim.range.umaps[[1]] <- ElbowPlot(Samples.combined,50)
for (dim in dim.range) {
  message(paste0("Testing step for dim: ",dim))
  Samples.combined <- RunUMAP(Samples.combined, reduction = "harmony", dims = 1:dim)
  dim.range.umaps[[paste0("dim:",dim)]]<-DimPlot(Samples.combined, reduction="umap", group.by = "Annotation.Level.3") + 
    NoLegend() + labs(title=paste0("dim: 1:",dim))+theme(axis.title.y = element_blank(),axis.text = element_text(size = 8),axis.title = element_text(size = 10))
}
print(ggarrange(plotlist = dim.range.umaps, widths = 20, heights = 30,ncol = 3,nrow = 3))
dev.off()
# Choose Dimension 1:30
Samples.combined <- RunUMAP(Samples.combined, reduction="harmony", dims=1:30, reduction.name = "MPCUMAP")
Samples.combined <- FindNeighbors(Samples.combined, reduction = 'harmony', dims = 1:30, verbose = FALSE)

#Step 4: clustering 
#Clustertree
sc.int=Samples.combined
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_02_cluster_resolution_tree.jpeg"), width=10 , height = 10)
#clusters and featureDotPlot for multiple resolution
message(paste0("start tryout for resolutions sample: ",sample))
pdf(file = paste0(sample, '_03_resolution tryout.pdf'),paper = "a4",height = 23,width = 33)
message(paste0("start loop for res sample: ",sample))
res.range = seq(0.2,0.6,by=0.2)
sc.int = Samples.combined
for (res in res.range) {
  sc.int <- FindClusters(sc.int, resolution = res) #cluster with res
  DefaultAssay(sc.int) <- "RNA"
  #plot with res on umap
  p1 <- DimPlot(sc.int, reduction = "MPCUMAP", label = TRUE, pt.size = 1, label.size = 6) + 
    NoLegend() + labs(title=paste0("res: ",res))+theme(axis.title.y = element_blank(),axis.text = element_text(size = 8),axis.title = element_text(size = 10))
  #plot nfeatures and ncount per cluster to identify low quality clusters
  v2 <- VlnPlot(sc.int,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0.1) +NoLegend() +theme(axis.text.x = element_text(angle = 0,hjust = 0),axis.text = element_text(size = 8),axis.title = element_text(size = 8),plot.title = element_text(size = 10))
  v3 <- VlnPlot(sc.int,group.by = "seurat_clusters", features =c('nCount_RNA'), pt.size = 0.1) +NoLegend() +theme(axis.text.x = element_text(angle = 0,hjust = 0),axis.text = element_text(size = 8),axis.title = element_text(size = 8),plot.title = element_text(size = 10))
  #top10 genes per cluster for res
  all.markers = FindAllMarkers(sc.int, test.use = "MAST",min.pct = 0.3,assay = "RNA")
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene)
  top10 <-unique(top10)
  dp1 <- DotPlot(sc.int, features = top10,dot.scale = 5,assay = "RNA")+       
    coord_flip() + NoLegend()+ theme(axis.title = element_blank(),axis.text.y=element_text(size = 7))
  print(ggarrange(plotlist = list(ggarrange(p1,v2,v3,ncol=1,nrow = 3),dp1),ncol = 2,widths = c(1.5,2)))
}
remove(sc.int)
dev.off()
message(paste0("DONE with sample: ",sample))
# Select res of 0.4 and 0.6 and annotate clusters based on marker expression
Samples.combined <- FindClusters(Samples.combined, resolution = 0.6) 
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
new.cluster.ids <- c("CD14 Mono #1", "FCGR3A Mono", "Res Mac", "CD14 Monocytes #2", "CD14 Monocytes #3","cDC #1", "SPP1 Mac","cDC #2", "T- & NK-Cells","EC", "cDC #3", "Fibro", "Prolif. cDC", "Prolif. MPC")
Samples.combined$Annotation_Level_2<- plyr::mapvalues(x=Samples.combined$RNA_snn_res.0.6, from = current.cluster.ids, to=new.cluster.ids)
Samples.combined <- FindClusters(Samples.combined, resolution = 0.4) 
current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
new.cluster.ids = c("CD14 Mono #1", "FCGR3A Mono", "CD14 Mono #2","Res Mac", "cDC", "SPP1 MAC", "T- & NK-Cells", "EC", "Prolif. cDC", "CXCL9 cDC", "Fibro", "Prolif. MPC") 
Samples.combined$Annotation_Level_1<- plyr::mapvalues(x=Samples.combined$RNA_snn_res.0.4, from = current.cluster.ids, to=new.cluster.ids)

#Remove all non MPC and proliferating clusters
Idents(Samples.combined) <- "Annotation_Level_2"
Samples.combined <- subset(Samples.combined, idents=c("T- & NK-Cells","EC","Fibro", "Prolif. cDC", "Prolif. MPC"), invert=TRUE)
Idents(Samples.combined) <- "Annotation_Level_1"
Samples.combined <- subset(Samples.combined, idents=c("T- & NK-Cells","Prolif. cDC"), invert=TRUE)
Samples.combined$Annotation_Level_1<- as.factor(Samples.combined$Annotation_Level_1)
Samples.combined$Annotation_Level_2<- as.factor(Samples.combined$Annotation_Level_2)
saveRDS(Samples.combined, file = paste0(sample, "05_seurat_object.RDS"))

# Dimplot of MPC
DimPlot(Samples.combined, reduction = "MPCUMAP", pt.size = 0.8)
ggsave(filename = paste0(sample,"04_Dimplot.jpeg"), width=6 , height = 5)

#Count number of cells
After_Subset<-length(Samples.combined@meta.data$orig.ident)
write.csv(After_Subset, file=paste0(sample,"_Cellnumber.csv"), row.names = F)

# Compositional Analysis with log2FC plots
Idents(Samples.combined)<- "Annotation_Level_1"
identities <- levels(Samples.combined@active.ident)
my_color_palette <- hue_pal()(length(identities))
dittoBarPlot(Samples.combined, var = "Annotation_Level_1", group.by = "Kidney.Function", scale = "percent")
dittobar_numbers=dittoBarPlot(Samples.combined, var="Annotation_Level_1", group.by="Kidney.Function", data.out = TRUE)
Cell_L1<-as.data.frame(dittobar_numbers$data)
Cell_L1_Healthy<- filter(Cell_L1, grouping=="Healthy")      %>% as.data.frame() %>% select(label, percent) %>% rename(Healthy="percent")
Cell_L1_CKD<-  filter(Cell_L1, grouping=="CKD")       %>% as.data.frame() %>% select(label, percent) %>% rename(CKD="percent")
Cell_L1 = merge(Cell_L1_Healthy,Cell_L1_CKD, by = "label")
Cell_L1$Log2FC <- log2(Cell_L1$CKD/Cell_L1$Healthy)
Cell_L1$label <- factor(Cell_L1$label, levels = Cell_L1$label[order(Cell_L1$Log2FC)])
ggplot(Cell_L1, aes(x=label, y=Log2FC)) + geom_bar(stat = "identity", fill=my_color_palette) + coord_flip()
ggsave(filename=paste0(sample, "06_WT_Cell_Log2FC.jpeg"), width= 3, height=5) + theme(text = element_text(size = 20,family = "Arial"), axis.text = element_text(size=20, family = "Arial"))
ggsave(filename=paste0(sample, "06_WT_Cell_Log2FC.svg"), width= 3, height=6)  + theme(text = element_text(size = 20,family = "Arial"), axis.text = element_text(size=20, family = "Arial"))

#Top3 Marker Genes as Dotplot
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
write_xlsx(all.markers, path=paste0(sample, "07_allmarkers_Level1.xlsx"), col_names = TRUE)
top3 <- all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top3 <-unique(top3)
DotPlot(Samples.combined, features = top3,dot.scale = 3, group.by = "Annotation_Level_1") +
  theme(axis.title = element_blank(), 
  text = element_text(size=8, family = "Arial"),
   axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
   axis.text.y = element_text(size=8, family = "Arial"),
  legend.text = element_text(size=6, family = "Arial"))
ggsave(filename = paste0(sample,"08_top3_featureDotplot_Level1_t.jpeg"), bg="white", width=5 , height = 2)
ggsave(filename = paste0(sample,"08_top3_featureDotplot_Level1_t.svg"), width=5 , height = 2)

#Feature Plots for selected markers
features_all = c("SPP1", "APOE")
Idents(Samples.combined)<- "Annotation_Level_1"
for (f in features_all) {
      FeaturePlot(Samples.combined, features =f, pt.size = 1, order=TRUE, label.size = 3.5,repel = TRUE, reduction = "MPCUMAP")
      ggsave(filename = paste0(sample,"_03_",f,".jpeg"), width=7 , height = 7)
      VlnPlot(Samples.combined, features = f,  group.by="Annotation_Level_1", flip = TRUE, pt.size = 0) +NoLegend() + theme(
        axis.title = element_blank(),
        text = element_text(size=8, family = "Arial"),
        axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
        axis.text.y = element_text(size=6, family = "Arial"),
        plot.title = element_text(size=9, family="Arial", face = "plain"),
        axis.line = element_line(size = 0.2), 
        axis.ticks = element_line(size = 0.2))
      ggsave(filename = paste0(sample,"_03_",f, "_Res02_Vln.jpeg"), width=2 , height = 2.2)
      ggsave(filename = paste0(sample,"_03_",f, "_Res02_Vln.svg"), width=2 , height = 2.2)
      }
