library(Seurat)
library(zellkonverter)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(tibble)
library(dittoSeq)
library(writexl)
# set fonts for ggplot
windowsFonts("Arial" = windowsFont("Arial"))

setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/integrated/")
sample="Cxcl4_snRNA_scvi_"

# Step 1: Import scvi integrated dataset and convert to seurat object
Integrated.5HAD <- readH5AD("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/snRNA-Data/for_sikander_konrad/scvi_output_clustered.h5ad")
Samples.combined <- as.Seurat(Integrated.5HAD, counts = "X", data = "X")
rm(Integrated.5HAD)

# Step 2: Normalisation, VST, Scaling
Samples.combined <- NormalizeData(Samples.combined, normalization.method = "LogNormalize", scale.factor = 10000)
Samples.combined <- FindVariableFeatures(Samples.combined, selection.method = "vst", nfeatures = 2000)
Samples.combined = ScaleData(Samples.combined, verbose = FALSE, features = rownames(Samples.combined))

# Adjust condition names
new.condition.names <- c("Cxcl4-/- Sham", "Cxcl4-/- IRI", "WT Sham", "WT IRI")
old.contidtion.names <- levels(Samples.combined$condition)
Samples.combined$condition<- plyr::mapvalues(x=Samples.combined$condition, from = old.contidtion.names, to=new.condition.names)
Samples.combined$condition<-factor(Samples.combined$condition, levels=c("WT Sham", "Cxcl4-/- Sham", "WT IRI", "Cxcl4-/- IRI"))
rm(new.condition.names, old.contidtion.names)

# Step 3: Check  annotated cell clusters
DimPlot(Samples.combined, group.by = "leiden0.5", reduction = "X_scvi_umap")
#Remove Adipocytes as these are artifacs from the capsule from isolation
Idents(Samples.combined)<-"leiden0.5"
Samples.combined<-subset(Samples.combined, idents = "10: Adipocytes", invert=TRUE)
Samples.combined$leiden0.5<-factor(Samples.combined$leiden0.5)
# Correct clusternames
new.leiden0.5.names <- c("DCT", "PT", "TAL", "Endo", "Fibro 1", "IC", "PC", "Leuko", "DTL", "Podo", "Injured Tub.", "VSMC", "PC/DME", "Neurons")
old.leiden0.5.names <- levels(Samples.combined$leiden0.5)
Samples.combined$leiden0.5<- plyr::mapvalues(x=Samples.combined$leiden0.5, from = old.leiden0.5.names, to=new.leiden0.5.names)
rm(new.leiden0.5.names, old.leiden0.5.names)

#Count number of cells
After_Subset<-length(Samples.combined@meta.data$orig.ident)
write.csv(After_Subset, file=paste0(sample,"_Cellnumber.csv"), row.names = F)

# Step 4: Dimplot and log2FC Compositional Analysis
# Dimplot
DimPlot(Samples.combined, group.by = "leiden0.5", reduction = "X_scvi_umap", pt.size = 0.3,split.by = "condition") + NoLegend()
ggsave(filename=paste0(sample, "02_Dimplot_NoLegend_Split.jpg"), width = 15, height =  6)
DimPlot(Samples.combined, group.by = "leiden0.5", reduction = "X_scvi_umap", label = TRUE, label.size = 4)
ggsave(filename=paste0(sample, "02_Dimplot.jpg"), width = 5, height =  5)

# Compositional Log2FC Plots
identities <- levels(Samples.combined@active.ident)
my_color_palette <- hue_pal()(length(identities))
dittoBarPlot(Samples.combined, "leiden0.5", group.by="condition", scale="percent", data.out=FALSE)
dittobar_numbers=dittoBarPlot(Samples.combined, "leiden0.5", group.by="condition", scale="percent", data.out = TRUE)
Cell_L1<-as.data.frame(dittobar_numbers$data)
Cell_L1_WTSHAM<- filter(Cell_L1, grouping=="WT Sham")      %>% as.data.frame() %>% select(label, percent) %>% rename(WT_SHAM="percent")
Cell_L1_WTIRI<-  filter(Cell_L1, grouping=="WT IRI")       %>% as.data.frame() %>% select(label, percent) %>% rename(WT_IRI="percent")
Cell_L1_KOSHAM<- filter(Cell_L1, grouping=="Cxcl4-/- Sham") %>% as.data.frame() %>% select(label, percent) %>% rename(Cxcl4KO_SHAM="percent")
Cell_L1_KOIRI<-  filter(Cell_L1, grouping=="Cxcl4-/- IRI")  %>% as.data.frame() %>% select(label, percent) %>% rename(Cxcl4KO_IRI="percent")
Cell_L1 = merge(Cell_L1_WTSHAM,Cell_L1_WTIRI, by = "label") %>% merge(Cell_L1_KOSHAM,by = "label")%>% merge(Cell_L1_KOIRI,by = "label")
Cell_L1$WT_Log2FC <- log2(Cell_L1$WT_IRI/Cell_L1$WT_SHAM)
Cell_L1$KO_Log2FC <- log2(Cell_L1$Cxcl4KO_IRI/Cell_L1$Cxcl4KO_SHAM)
Cell_L1$label <- factor(Cell_L1$label, levels = Cell_L1$label[order(Cell_L1$WT_Log2FC)])
ggplot(Cell_L1, aes(x=label, y=WT_Log2FC)) + geom_bar(stat = "identity", fill=my_color_palette) + coord_flip() + ylim(-2,3) +
  theme(text = element_text(size=8, family = "Arial"),
        axis.text.x = element_text(size=6, family = "Arial"),
        axis.text.y = element_text(size=8, family = "Arial"),
        axis.title = element_text(size=8, family = "Arial"),
        plot.title = element_text(size=9, family="Arial", face = "plain"))
ggsave(filename=paste0(sample, "03_WT_Cell_Log2FC.jpeg"), width= 1.5, height=2)
ggsave(filename=paste0(sample, "03_WT_Cell_Log2FC.svg"), width= 1.5, height=2)
ggplot(Cell_L1, aes(x=label, y=KO_Log2FC)) + geom_bar(stat = "identity", fill=my_color_palette) + coord_flip() + ylim(-2,3) +scale_x_discrete(position = "top")+
  theme(text = element_text(size=8, family = "Arial"),
        axis.text.x = element_text(size=6, family = "Arial"),
        axis.text.y = element_text(size=8, family = "Arial"),
        axis.title = element_text(size=8, family = "Arial"),
        plot.title = element_text(size=9, family="Arial", face = "plain"))
ggsave(filename=paste0(sample, "04_KO_Cell_Log2FC.jpeg"), width= 1.5, height=2)
ggsave(filename=paste0(sample, "04_KO_Cell_Log2FC.svg"), width= 1.5, height=2)

# Step 5: Top3 Upregulated Genes as Dotplot
Idents(Samples.combined)<-"leiden0.5"
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "originalexp")
write_xlsx(all.markers, path=paste0(sample, "05_allmarkers_Level1.xlsx"), col_names = TRUE)
top3 <- all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top3 <-unique(top3)
DotPlot(Samples.combined, features = top3,dot.scale = 4, group.by = "leiden0.5") + coord_flip() + theme(
  axis.title = element_blank(),
  text = element_text(size=8, family = "Arial"),
  axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=8, family = "Arial"))
ggsave(filename = paste0(sample,"06_top3_featureDotplot.jpeg"), width=4 , height = 7)
ggsave(filename = paste0(sample,"06_top3_featureDotplot.svg"), width=4 , height = 7)

# Step 6: Featureplots of interesting genes
features=c("Spp1", "Apoe", "C1qa", "Arg1", "Myh11", "CD74", "Ptprc")
for (f in features) {
  FeaturePlot(Samples.combined, features = f, min.cutoff = "q9", order=TRUE, pt.size = 1, reduction = "X_scvi_umap") 
  ggsave(filename = paste("FeatureandVlnPlots/", sample,"Feature_",f,".jpeg"), width=10 , height = 7)
  VlnPlot(Samples.combined, features = f, group.by="leiden0.5", pt.size = 0)+ NoLegend()+ theme(
    axis.title = element_blank(),
    text = element_text(size=8, family = "Arial"),
    axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
    axis.text.y = element_text(size=8, family = "Arial"),
    plot.title = element_text(size=9, family="Arial", face = "plain"),
    axis.line = element_line(size = 0.2), 
    axis.ticks = element_line(size = 0.2))
  ggsave(filename = paste("FeatureandVlnPlots/",sample,"Vln_",f,".jpeg"), width=3.3 , height = 2)
  ggsave(filename = paste("FeatureandVlnPlots/",sample,"Vln_",f,".svg"), width=3.3 , height = 2)
  }
