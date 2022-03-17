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
library(tibble)
library(pheatmap)
library(readxl)

#Working directory
# Load MPC Sample
setwd("c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/Forte_MI_data/ImmuneCellAnalysis/")
sample="MI_ImmuneCell"
Samples.combined <- readRDS(file = "c:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/Forte_MI_data/ImmuneCellAnalysis/MI_ImmuneCell_01_SeuratObj_AllImmuneCells.RDS") #if its an rds file
source(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/Volcano_plot_KH.R")
source_python("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/huCKD_Nature_Decode_Myofibroblast_Dataset/HuPDGFRb/Correlation Analysis/Python/Python-Correlation.py")

pdf(file = paste0(sample, '_02_dimension tryout.pdf'),paper = "a4",height = 23,width = 33)
dim.range = seq(5,50,by=5)
dim.range.umaps = vector("list")
dim.range.umaps[[1]] <- ElbowPlot(Samples.combined,50)
for (dim in dim.range) {
  message(paste0("Testing step for dim: ",dim))
  Samples.combined <- RunUMAP(Samples.combined, dims = 1:dim)
  dim.range.umaps[[paste0("dim:",dim)]]<-DimPlot(Samples.combined, group.by = "seurat_clusters") + 
    NoLegend() + labs(title=paste0("dim: 1:",dim))+theme(axis.title.y = element_blank(),axis.text = element_text(size = 8),axis.title = element_text(size = 10))
}
print(ggarrange(plotlist = dim.range.umaps, widths = 20, heights = 30,ncol = 3,nrow = 3))
dev.off()

#Choose Dimensions 1:25 based on data
DefaultAssay(Samples.combined ) <- "integrated"
Samples.combined <- RunUMAP(Samples.combined, dims=1:25)
Samples.combined <- FindNeighbors(Samples.combined, reduction = 'pca', dims = 1:25, verbose = TRUE)

# Resolution tryout
message(paste0("start tryout for resolutions sample: ",sample))
pdf(file = paste0(sample, '_05_resolution tryout.pdf'),paper = "a4",height = 23,width = 33)
#clusters and featureDotPlot for multiple resolution
message(paste0("start loop for res sample: ",sample))
res.range = seq(0.2,1.3,by=0.1)
for (res in res.range) {
  DefaultAssay(Samples.combined) <- "integrated"
  Samples.combined <- FindClusters(Samples.combined, resolution = res) #cluster with res
  DefaultAssay(Samples.combined) <- "RNA"
  #plot with res on umap
  p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) + 
    NoLegend() + labs(title=paste0("res: ",res))+theme(axis.title.y = element_blank(),axis.text = element_text(size = 8),axis.title = element_text(size = 10))
  #plot nfeatures and ncount per cluster to identify low quality clusters
  v2 <- VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0.1) +NoLegend() +theme(axis.text.x = element_text(angle = 0,hjust = 0),axis.text = element_text(size = 8),axis.title = element_text(size = 8),plot.title = element_text(size = 10))
  v3 <- VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nCount_RNA'), pt.size = 0.1) +NoLegend() +theme(axis.text.x = element_text(angle = 0,hjust = 0),axis.text = element_text(size = 8),axis.title = element_text(size = 8),plot.title = element_text(size = 10))
  #top10 genes per cluster for res
  all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene)
  top10 <-unique(top10)
  dp1 <- DotPlot(Samples.combined, features = top10,dot.scale = 5,assay = "RNA")+       
    coord_flip() + NoLegend()+ theme(axis.title = element_blank(),axis.text.y=element_text(size = 7))
  print(ggarrange(plotlist = list(ggarrange(p1,v2,v3,ncol=1,nrow = 3),dp1),ncol = 2,widths = c(1.5,2)))
}
dev.off()
message(paste0("DONE with sample: ",sample))

# Select Res 0.2 and 0.9 for Annotation
Samples.combined <- FindClusters(Samples.combined, resolution = 0.2)
current.cluster.ids = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
new.cluster.ids = c("Res Mac","Granulocytes", "Spp1 Mac","B-Cells","Ly6c2hi Mono","cDC2", "T- & NK-Cells", "Ifn Mac", "Fibroblasts","B-Cell/EC Doublets")
Samples.combined$Annotation_Level_1<- plyr::mapvalues(x=Samples.combined$integrated_snn_res.0.2, from = current.cluster.ids, to=new.cluster.ids)
Samples.combined <- FindClusters(Samples.combined, resolution = 0.9)
current.cluster.ids = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21")
new.cluster.ids = c("Ly6c2hi Mono", "B-Cells #1", "Granulocytes #1", "Res Mac #1", "Res Mac #2", "Res Mac #3", "Res Mac #4", "cDC2", "Granulocytes #2", "Mrc1 Mac", "T-Cells", "Ifn Mac", "Ccr2 Mac #2", "Spp1 Mac #1", "B Cells #2", "Spp1 Mac #2", "NK-Cells", "Fibroblast/MPC Doublets", "Prol. cDC2", "Ccr7 DC", "B-Cell/EC Doublets", "Granulocyte/EC Doublets")
Samples.combined$Annotation_Level_2<- plyr::mapvalues(x=Samples.combined$seurat_clusters, from = current.cluster.ids, to=new.cluster.ids)

# Remove non Immune Cell Clusters
Samples.combined = subset(Samples.combined, idents = c(0,1,2,3,4,5,6,7,9,8,10,11,12,13,14,15,16,19))

#Add timepoint and remove timepoint "None" as identity is unclear
Idents(Samples.combined)<-"orig.ident"
current.cluster.ids = levels(Idents(object = Samples.combined))
new.cluster.ids = c("Sham", "Sham", "None", "MI_Day1", "MI_Day3", "MI_Day5", "MI_Day7", "MI_Day14", "MI_Day28")
Samples.combined$timepoint<- plyr::mapvalues(x=Samples.combined$orig.ident, from = current.cluster.ids, to=new.cluster.ids)
Idents(Samples.combined)<-"timepoint"
Samples.combined = subset(Samples.combined, idents = c("Sham", "MI_Day1", "MI_Day3", "MI_Day5", "MI_Day7", "MI_Day14", "MI_Day28"))

After_QC<-length(Samples.combined@meta.data$orig.ident)
write.csv(After_QC, file=paste0(sample,"_Cellnumber.csv"), row.names = F)

DefaultAssay(Samples.combined ) <- "integrated"
Samples.combined <- RunUMAP(Samples.combined, dims=1:25)
Samples.combined <- FindNeighbors(Samples.combined, reduction = 'pca', dims = 1:25, verbose = TRUE)
Idents(Samples.combined)<-"Annotation_Level_1"

# Dimplot with Annotation
DimPlot(Samples.combined, reduction = "umap", group.by = "Annotation_Level_1" , pt.size = 0.8)
ggsave(filename=paste0(sample, "_06_Level1_Res0.2.jpeg"),width=9, height= 8)

## Feature and Count plots
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_08_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "Annotation_Level_2", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_09_nfeature_ncount_vln.jpeg"), width=10 , height =5)

#Top3 Marker Genes as Dotplot
all.markers <- FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
write_xlsx(all.markers, path=paste0(sample, "_10_allmarkers_Level1.xlsx"), col_names = TRUE)
all.markers <- read_excel(path = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/Forte_MI_data/ImmuneCellAnalysis/MI_ImmuneCell_10_allmarkers_Level1.xlsx")
top3 <- all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top3 <-unique(top3)
DotPlot(Samples.combined, features = top3,dot.scale = 3, group.by = "Annotation_Level_1") + theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(
  axis.title = element_blank(),
  text = element_text(size=8, family = "Arial"),
  axis.text.x = element_text(size=7, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=8, family = "Arial"))
ggsave(filename = paste0(sample,"_10_top3_featureDOTplot_Level1.jpeg"), bg="white", width=5.1, height = 2)
ggsave(filename = paste0(sample,"_10_top3_featureDOTplot_Level1.svg"), width=5.1, height = 2)

# Compositional Analysis with Dittobarplot
dittoBarPlot(Samples.combined, "Annotation_Level_1", group.by="timepoint", scale="percent", data.out=FALSE, x.reorder = c(7,1,4,5,6,2,3))
ggsave(filename = paste0(sample,"_12_Dittobarplot.jpeg"), width=6, height = 5)
ggsave(filename = paste0(sample,"_12_Dittobarplot.svg"), width=5, height = 3)

#Save Seurat Object
saveRDS(Samples.combined, file = paste0(sample, "_13_SeuratObjt_AfterFirstAnalysis.RDS"))

# Calculate ECM Regulator and Core Matrisome Scores
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

DefaultAssay(Samples.combined)="RNA"
ctrl_genes = 50 #important
gsets=c("ECM_Regulators","Core_matrisome")

for (gset in gsets){
  features = matrisome_mm_genesetlist[gset]
  message(gset)
  # Add average expression of genes in gset minus the average expression of
  # 35 random genes 
  Samples.combined = AddModuleScore(object = Samples.combined, features = features, name = gset, ctrl = ctrl_genes)
  FeaturePlot(Samples.combined, features = paste0(gset, '1'),pt.size = 1,order = T,label = FALSE) + 
          scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red')
  ggsave(filename = paste0(sample,"_",gset,"_f.jpeg"), width=10 , height = 10)
  ggsave(filename = paste0(sample,"_",gset,"_f.svg"), width=10 , height = 10)
  VlnPlot(Samples.combined, features = paste0(gset, '1'),  group.by="Annotation_Level_1",pt.size = 0, sort = "decreasing")
    +NoLegend() + coord_flip() +geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)
  ggsave(filename = paste(sample,"_", gset, "_V.jpeg"), width=4, height = 5.5)
  ggsave(filename = paste(sample,"_", gset, "_V.svg"), width=4, height = 5.5)
}

#ttest_calculations for ECM Regulator Score
Naba.df<-FetchData(object = Samples.combined, vars=c("Annotation_Level_1", "ECM_Regulators1", "Pf4"))
sink(paste0("",sample,"_14_TTest_Results.txt"))
print(paste0("t-test for  ECM Regulator Spp1 Mac WT IRI vs Res Mac"))
print(t.test(Naba.df[which(Naba.df$Annotation_Level_1 =="Spp1 Mac"),2], Naba.df[which(Naba.df$Annotation_Level_1=="Res Mac"),2], alternative = "two.sided"))
sink()

# Perform gene correlation to ECM Regulator genes to find Top Correlating Genes with ECM-Regulator Score
Samples.combined@assays$ECMRegulators <- CreateAssayObject(counts=as.data.frame(Samples.combined$ECM_Regulators1))
# Cluster at a high resolution for later Pseudobulking for correlation to reduce data sparsity for correlation analysis
Samples.combined <- FindClusters(Samples.combined, resolution = 5)

# Remove ECM Regulating Genes for Correlation Analysis
ECMRegulator_Genes <- matrisome_mm_genesetlist$ECM_Regulators
counts <- GetAssayData(Samples.combined, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% ECMRegulator_Genes)),]
Samples.combined_minusECMReg <- subset(Samples.combined, features = rownames(counts))
MeanExpr<-as.data.frame(t(as.data.frame(AverageExpression(Samples.combined_minusECMReg, group.by = "seurat_clusters", assays = "RNA"))))

# Correlation
MeanExprECMReg = as.data.frame(t(as.data.frame(AverageExpression(Samples.combined, group.by = "seurat_clusters", assays = "ECMRegulators"))))
colnames(MeanExprECMReg)<-"ECMRegulators"
CorMatrix=cbind(MeanExprECMReg, MeanExpr)
cors=as.data.frame(py$correlate_means_to_gene(CorMatrix, corr_gene = "ECMRegulators"))
color_heatmap = colorRampPalette(c("darkblue", "whitesmoke", "indianred"))(100)
cors<-rownames_to_column(cors, var="gene")
Top40_Up<-head(cors[,1], n=21L)
MeanExpr_Sorted=as.data.frame(t(as.data.frame(AverageExpression(Samples.combined, group.by = "seurat_clusters", assays="RNA",features = Top40_Up))))
MeanExpr_Sorted=as.data.frame(cbind(MeanExprECMReg, MeanExpr_Sorted)) 
MeanExpr_Sorted <- MeanExpr_Sorted%>% arrange((ECMRegulators)) %>% t
pheatmap(MeanExpr_Sorted, color = color_heatmap, cluster_cols = FALSE, cluster_rows = FALSE,
           scale = "row", show_colnames = FALSE, 
           filename = paste0(sample,"_ECMRegulators_Top20_Cor.pdf"), width= 3, height = 3)

# Plot Cxcl4 Expression
print(FeaturePlot(Samples.combined, features = "Pf4",pt.size = 1,order = T,label = FALSE) + 
        scale_colour_gradient2(low = 'blue', mid =  'lightgrey', high = 'red'))
ggsave(filename = paste0("MI_dataAll_03_FeatureAndVlnPlots/",sample,"_Pf4_F.jpeg"), width=10 , height = 10)
VlnPlot(Samples.combined, features = "Pf4",  group.by="Annotation_Level_1",pt.size = 0, sort = "decreasing")+NoLegend() + coord_flip() +geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)
ggsave(filename = paste("MI_dataAll_03_FeatureAndVlnPlots/", sample,"_Pf4_V.jpeg"), width=4, height = 5.5)
ggsave(filename = paste("MI_dataAll_03_FeatureAndVlnPlots/", sample,"_Pf4_V.svg"), width=4, height = 5.5)

# Plot Trem2 Expression
FeaturePlot(Samples.combined, features = "Ccr2",pt.size = 1,order = T,label = FALSE) + scale_colour_gradient2(low = 'blue', mid =  'lightgrey', high = 'red')
ggsave(filename = paste0("MI_dataAll_03_FeatureAndVlnPlots/",sample,"_Pf4_F.jpeg"), width=10 , height = 10)
VlnPlot(Samples.combined, features = "Ccr2",  group.by="Annotation_Level_1",pt.size = 0, sort = "decreasing")+NoLegend() + coord_flip() +geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)
ggsave(filename = paste("MI_dataAll_03_FeatureAndVlnPlots/", sample,"_Pf4_V.jpeg"), width=4, height = 5.5)
ggsave(filename = paste("MI_dataAll_03_FeatureAndVlnPlots/", sample,"_Pf4_V.svg"), width=4, height = 5.5)

# Violin Plots for Marker genes Fn1, Spp1 and C1qa for Figure S2E
features=c("Fn1", "Spp1", "C1qa")
for (f in features) {
  VlnPlot(Samples.combined, features = f,  group.by="Annotation_Level_1", pt.size = 0) + NoLegend()+ 
    geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)+theme(
      axis.title = element_blank(),
      text = element_text(size=8, family = "Arial"),
      axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
      axis.text.y = element_text(size=8, family = "Arial"),
      plot.title = element_text(size=9, family="Arial", face = "plain"),
      axis.line = element_line(size = 0.2), 
      axis.ticks = element_line(size = 0.2))
  ggsave(filename = paste(sample,"Vln_",f,".jpeg"), width=3, height = 1.7)
  ggsave(filename = paste(sample,"Vln_",f,".svg"), width=3, height = 1.7)
}