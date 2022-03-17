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

setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/AllIntegrated/")

sample = 'Rao'

#Step 1: Read in All Data
# Read in all data 
# Of note: the data-sets ICM2 and ICM3 were mislabeled as RVP (=NMI) or LVP (=MI) on the GEO database. We corrected this after reading in the data.
# QC metrics were kept the same as in Rao et al.
N1_LVP <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/N-1/N1_LVP/")
N1_LVP <- CreateSeuratObject(counts = N1_LVP, project = "N1_LVP", min.cells = 3, min.features = 500)
N1_LVP$patient.id <- "N1"
N1_LVP$disease <- "Healthy"
N1_LVP$location <- "LVP"

N1_RVP <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/N-1/N1_RVP/")
N1_RVP <- CreateSeuratObject(counts = N1_RVP, project = "N1_RVP", min.cells = 3, min.features = 500)
N1_RVP$patient.id <- "N1"
N1_RVP$disease <- "Healthy"
N1_RVP$location <- "RVP"

ICM1_MI <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/ICM-1/ICM1_MIP/")
ICM1_MI <- CreateSeuratObject(counts = ICM1_MI, project = "ICM1_MI", min.cells = 3, min.features = 500)
ICM1_MI$patient.id <- "ICM1"
ICM1_MI$disease <- "ICM"
ICM1_MI$location <- "MI"

ICM1_NMI <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/ICM-1/ICM1_NMIP/")
ICM1_NMI <- CreateSeuratObject(counts = ICM1_NMI, project = "ICM1_NMI", min.cells = 3, min.features = 500)
ICM1_NMI$patient.id <- "ICM1"
ICM1_NMI$disease <- "ICM"
ICM1_NMI$location <- "NMI"

ICM2_MI <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/ICM-2/ICM2_LVP/")
ICM2_MI <- CreateSeuratObject(counts = ICM2_MI, project = "ICM2_MI ", min.cells = 3, min.features = 500)
ICM2_MI$patient.id <- "ICM2"
ICM2_MI$disease <- "ICM"
ICM2_MI$location <- "MI"

ICM2_NMI <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/ICM-2/ICM2_RVP/")
ICM2_NMI <- CreateSeuratObject(counts = ICM2_NMI, project = "ICM2_NMI ", min.cells = 3, min.features = 500)
ICM2_NMI$patient.id <- "ICM2"
ICM2_NMI$disease <- "ICM"
ICM2_NMI$location <- "NMI"

ICM3_MI <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/ICM-3/ICM3_LVP/")
ICM3_MI <- CreateSeuratObject(counts = ICM3_MI, project = "ICM3_MI", min.cells = 3, min.features = 500)
ICM3_MI$patient.id <- "ICM3"
ICM3_MI$disease <- "ICM"
ICM3_MI$location <- "MI"

ICM3_NMI <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/ICM-3/ICM3_RVP/")
ICM3_NMI <- CreateSeuratObject(counts = ICM3_NMI, project = "ICM3_NMI ", min.cells = 3, min.features = 500)
ICM3_NMI$patient.id <- "ICM3"
ICM3_NMI$disease <- "ICM"
ICM3_NMI$location <- "NMI"

DCM2_LVP <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/DCM-2/DCM2_LVP/")
DCM2_LVP <- CreateSeuratObject(counts = DCM2_LVP, project = "DCM2_LVP", min.cells = 3, min.features = 500)
DCM2_LVP$patient.id <- "DCM2"
DCM2_LVP$disease <- "DCM"
DCM2_LVP$location <- "LVP"

DCM2_RVP <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/DCM-2/DCM2_RVP/")
DCM2_RVP <- CreateSeuratObject(counts = DCM2_RVP, project = "DCM2_RVP", min.cells = 3, min.features = 500)
DCM2_RVP$patient.id <- "DCM2"
DCM2_RVP$disease <- "DCM"
DCM2_RVP$location <- "RVP"

DCM3_LVP <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/DCM-3/DCM3_LVP/")
DCM3_LVP <- CreateSeuratObject(counts = DCM3_LVP, project = "DCM3_LVP", min.cells = 3, min.features = 500)
DCM3_LVP$patient.id <- "DCM3"
DCM3_LVP$disease <- "DCM"
DCM3_LVP$location <- "LVP"

DCM3_RVP <- Read10X(data.dir = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/SelectedDatasetsforAnalysis/DCM-3/DCM3_RVP/")
DCM3_RVP <- CreateSeuratObject(counts = DCM3_RVP, project = "DCM3_RVP", min.cells = 3, min.features = 500)
DCM3_RVP$patient.id <- "DCM3"
DCM3_RVP$disease <- "DCM"
DCM3_RVP$location <- "RVP"

#Merge Datasets
sc <- merge(N1_LVP, y=c(N1_RVP, ICM1_MI, ICM1_NMI, ICM2_MI, ICM2_NMI, ICM3_MI, ICM3_NMI, DCM2_LVP, DCM2_RVP, DCM3_LVP, DCM3_RVP), project="sc")
rm(N1_LVP, N1_RVP, ICM1_MI, ICM1_NMI, ICM2_MI, ICM2_NMI, ICM3_MI, ICM3_NMI, DCM2_LVP, DCM2_RVP, DCM3_LVP, DCM3_RVP)

# Step 2: QC_Check
#'^mt-' = genes from the features.tsv.gz file that start with mt-. This pattern should be lower case for mice and upper case for humans.
sc[['percent.mt']] = PercentageFeatureSet(sc, pattern = '^MT-')
VlnPlot(sc, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.1)
ggsave(filename = paste0(sample,"_01_Gene&RNA count Filtered.jpeg"), width=10 , height = 10)
mt_VlnPlot <- VlnPlot(sc, features = 'percent.mt', pt.size = 0.1) + geom_hline(yintercept = 10)
ggsave(filename = paste0(sample,"_01_mtVlnPlot with Cutoff.jpeg"), width=10 , height = 10)
nFeature_VlnPlot <- VlnPlot(sc, features = 'nCount_RNA', pt.size = 0.1) + geom_hline(yintercept = 8000) + geom_hline(yintercept = 800)
ggsave (filename = paste0(sample,"_01_nCount_RNA with Cutoff.jpeg"), width=10 , height = 10)

# Step 3: Quality filter - The Same Quality filter Metrics were used as by Rao et al. (min. features 500, nRNA <800 and >8000 excluded)
Before_QC = length(sc@meta.data$orig.ident)
sc = subset(sc, subset = nCount_RNA > 800 & nCount_RNA < 8000 & percent.mt <10)
After_QC = length(sc@meta.data$orig.ident)

# Normalize, find variable features & scale
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sc <- FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
sc <- ScaleData(sc, verbose=FALSE)

############################## -- HARMONY INTEGRATION --#######################################
# Dimplot before Integration
sc <- RunPCA(sc, ncps=30, verbose=FALSE)
sc <- RunUMAP(sc, reduction="pca", dims=1:30)
DimPlot(sc, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, ncol=4) + NoLegend()
ggsave(filename = paste0(sample, "_02_Dimplot_beforeInteg.jpeg"), width= 8, height = 8)

#Harmony-Integration
sc <- sc %>% RunHarmony("orig.ident", plot_convergence=TRUE)
ggsave(filename = paste0(sample, "_03_Harmony.jpeg"), width= 4, height =4)
ggsave(filename = paste0(sample, "_03_Harmony.pdf"), width= 4, height =4)

# Dimplot After Integration
sc <- RunUMAP(sc, reduction="harmony", dims=1:30)
sc <- FindNeighbors(sc, reduction = 'harmony', dims = 1:30, verbose = FALSE)
DimPlot(sc, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, ncol=2) + NoLegend()
ggsave(filename = paste0(sample, "_02_Dimplot_afterInteg.jpeg"), width= 8, height = 8)

# step 4: Save integrated dataset before clustering
saveRDS(sc, file = paste0(sample, "_04_seurat_obj_after_Int.rds"))
sc<-readRDS(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/in Human/Rao_MI_scRNA/AllIntegrated/Rao_04_seurat_obj_after_Int.rds")

# Step 5: Feature, Count and mt.gene Plots
FeaturePlot(sc, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_05_featureplot.jpeg"), width=20 , height = 10)

VlnPlot(sc,features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), pt.size = 0) & NoLegend() & theme(
  axis.title = element_blank(),
  text = element_text(size=8, family = "Arial"),
  axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
  axis.text.y = element_text(size=6, family = "Arial"),
  plot.title = element_text(size=9, family="Arial", face = "plain"),
  axis.line = element_line(size = 0.2), 
  axis.ticks = element_line(size = 0.2))
ggsave(filename = paste0(sample,"_05_Gene&RNA count Filtered.jpeg"), width=8 , height = 2)
ggsave(filename = paste0(sample,"_05_Gene&RNA count Filtered.svg"), width=8 , height = 2)

# Step 6: Find optimal dimension for UMAP
message(paste0("start tryout for resolutions sample: ",sample))
pdf(file = paste0(sample, '_06_dimension tryout.pdf'),paper = "a4",height = 23,width = 33)
#PCA testing
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

# choose dimension of 1:30
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:30)
sc <- FindNeighbors(sc, reduction = 'harmony', dims = 1:30, verbose = FALSE)

# Step 7: Clustering 
# Clustree
sc.int=sc
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_06_cluster_resolution_tree.jpeg"), width=10 , height = 10)

# Select resolution 0.1 for simple subsetting of MPC
sc <- FindClusters(sc, resolution=0.1)

DimPlot(sc, group.by="seurat_clusters",pt.size = 0.8)
ggsave(filename = paste0(sample, "_07_Dimplot.jpeg"), width=6 , height = 5)
DimPlot(sc, group.by="seurat_clusters",pt.size = 0.8, label = TRUE)
ggsave(filename = paste0(sample, "_07_Dimplot_l.jpeg"), width=6 , height = 5)

#Top3 Marker Genes as Dotplot at Level1
all.markers = FindAllMarkers(sc, test.use = "MAST",min.pct = 0.3,assay = "RNA")
write_xlsx(all.markers, path=paste0(sample, "07_allmarkers_Res01.xlsx"), col_names = TRUE)
top3 <- all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) %>% pull(gene)
top3 <-unique(top3)
DotPlot(sc, features = top3,dot.scale = 3, group.by = "seurat_clusters") +
  theme(axis.title = element_blank(), 
        text = element_text(size=8, family = "Arial"),
        axis.text.x = element_text(size=8, family = "Arial", angle=45, hjust=1),
        axis.text.y = element_text(size=8, family = "Arial"),
        legend.text = element_text(size=6, family = "Arial"))
ggsave(filename = paste0(sample,"_07_top3_featureDotplot_Level1_t.jpeg"), bg="white", width=5 , height = 2)
ggsave(filename = paste0(sample,"_07_top3_featureDotplot_Level1_t.svg"), width=5 , height = 2)

# Step 8: Subset MPC
# Subset MPC  for further analysis and save
sc.MPC <- subset(sc,idents = c(1,3))
MPC = length(sc.MPC@meta.data$orig.ident)
saveRDS(sc.MPC, file = paste0(sample, "_08_seurat_obj_MPC.rds"))

filter_step_counter=list()
filter_step_counter<-c(Before_QC, After_QC, MPC, Mac)
filter_step_counter<- as.data.frame(filter_step_counter)
row.names(filter_step_counter) = c("Before QC:", "After QC", "MPC", "Mac")
filter_step_counter<- as.data.frame(filter_step_counter)
filter_step_counter<-t(filter_step_counter)
filter_step_counter<- as.data.frame(filter_step_counter)
write_xlsx(x = filter_step_counter,path = paste0(sample,"_09_table of filter_step_counters.xlsx"))


