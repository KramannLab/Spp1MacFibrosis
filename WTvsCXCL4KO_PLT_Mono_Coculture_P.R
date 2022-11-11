# PROGGENY, DOROTHEA AND GO Analysis of BULK RNA Script
library(ggplot2)
library(reshape2)
library(writexl)
library(readxl)
library(plyr)
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)
library(readr)
library(ggrepel)
library(gridExtra)
library(graphics)
library(fgsea)
library(stringr)
library(Seurat)

setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/Manuscript/Final/Cell reports/Revision/BulkRNA/")
source(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/Pseudobulking_support_functions.R")
source(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/Volcano_plot_KH.R")
sc <- readRDS("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vivo/Forte_MI_data/ImmuneCellAnalysis/MI_ImmuneCell_15_SeuratObjt_AfterAnnotation.RDS")

DEG_KOvsWT<- as.data.frame(read_excel(path = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/2. CXCl4/In Vitro/PF4-KOvWT - mPBMC-PLT Coculture - BulkSeq 180821/Analysis/SupplementalTable7_WtvsCxcl4KO_PLTstim_Mono_DEG.xlsx", na = c(NA, "NA")))
colnames(DEG_KOvsWT)<-c("ID", "base.Meanx", "log2FC", "lfc.SE", "t", "pval", "adj.pval", "index")
df<- DEG_KOvsWT
sample <- "Plt_Mono_KOvsWT"


# Remove genes with NA for p.values as these are very weakly expressed and therefore assigned NA values by DESeq comparison
df<-na.omit(df)
volcano_nice(df = df, hAss = 0.5, FCIndex = 3,pValIndex = 6, IDIndex = 1, vAss = 0.1, label = TRUE, straight = FALSE, nlabels = 5, manual_labels = c("Arg1", "Spp1", "Fn1") ) +
  labs(x="logFC", y="-log(P. Val.)", title=sample) + theme(plot.title = element_text(hjust=0.5))
ggsave(filename = paste0(sample,"01_Vulcano.jpeg"), width=6 , height = 5, bg = "white", dpi = 600)
ggsave(filename = paste0(sample,"01_Vulcano.svg"), width=6 , height = 5, bg = "white", dpi = 600)

###### Progeny Analysis of TTOP-Results
paletteLength <-100
myColor <-colorRampPalette(c("darkblue", "whitesmoke", "indianred"))(paletteLength)

df_matrix <- df %>% 
  dplyr::select(ID,t) %>% 
  dplyr::filter(!is.na(t)) %>%
  column_to_rownames(var="ID")%>%  
  as.matrix()

PathwayActivity_counts = progeny(df_matrix, scale=FALSE, organism = "Mouse", top=500 , perm=10000,z_scores = TRUE)
rownames(PathwayActivity_counts) = sample
Activity_counts=as.vector(PathwayActivity_counts)
progenyBreaks = c(seq(min(Activity_counts),0,length.out=ceiling(paletteLength/2)+1), seq(max(Activity_counts)/paletteLength, max(Activity_counts), length.out = floor(paletteLength/2)))
PathwayActivity_zscore <- t(PathwayActivity_counts)
colnames(PathwayActivity_zscore)<-"NES"

PathwayActivity_zscore_df <- as.data.frame((PathwayActivity_zscore))%>%
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway=factor(Pathway))

ggplot(PathwayActivity_zscore_df,aes(y = reorder(Pathway, NES), x = NES)) + geom_bar(aes(fill = NES), stat = "identity") + labs(title = paste0(sample, " Progeny NES"))+
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal()+
  theme(plot.title =  element_text(face="bold", size = 14, hjust = 0.5), axis.title=element_text(face="bold", size=12),
        axis.text.x=element_text(angle=45, hjust=1, size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"),
        panel.grid.minor = element_blank()) +
  xlab("Pathways")
ggsave(filename = paste0(sample,"02_Progeny_NES.jpeg"), width=3, height=7, bg = "white", dpi = 600)
ggsave(filename = paste0(sample,"02_Progeny_NES.svg"), width=3, height=7, bg = "white")


#Infer Transcription Factor activity using Dorothea
data(dorothea_mm, package="dorothea")
regulons <- dorothea_mm %>% dplyr::filter(confidence %in% c("A", "B", "C"))

current_df_matrix <- df %>% 
    dplyr::select(ID,t) %>% 
    dplyr::filter(!is.na(t)) %>%
    column_to_rownames(var="ID")%>%
    as.matrix()

# CHECK WHETHER PLEIOTROPY WORKS OR NOT (#pleiotropy = TRUE)
tf_activities_stat = dorothea::run_viper(current_df_matrix, regulons, 
                                           options = list(minsize=30, eset.filter=FALSE,
                                                          cores=1, verbose=FALSE, nes=TRUE))

tf_activities_stat_top25<- tf_activities_stat%>%
    as.data.frame() %>%
    rownames_to_column(var="GeneID") %>%
    dplyr::rename(NES="t") %>%
    dplyr::top_n(20, wt=abs(NES)) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(GeneID=factor(GeneID))
  
ggplot(tf_activities_stat_top25, aes(y=reorder(GeneID, NES), x=NES)) + labs(title = paste0(sample, " Dorothea NES")) +
    geom_bar(aes(fill=NES), stat="identity") + 
    scale_fill_gradient2(low="darkblue", high="indianred", mid="whitesmoke", midpoint=0)+
    theme_minimal()+
    theme(plot.title =  element_text(face="bold", size = 14, hjust = 0.5), axis.title=element_text(face="bold", size=12),
          axis.text.x=element_text(angle=45, hjust=1, size=10, face="bold"),
          axis.text.y = element_text(size=10, face="bold"),
          panel.grid.minor = element_blank()) + 
    xlab("Transcription Factor")
ggsave(filename = paste0(sample, "03_Dorothea.jpeg"), width = 3, height = 6, bg = "white", dpi = 600)
ggsave(filename = paste0(sample, "03_Dorothea.svg"), width = 3, height = 6, bg = "white")
tf_xlsx<-as.data.frame(tf_activities_stat) %>% rownames_to_column()
colnames(tf_xlsx)<-c("TF", "NES")
write_xlsx(x = tf_xlsx,path = paste0(sample, "03_DorotheaTF.xlsx"))

#Score Forte Immune cells based on Top DE Genes
df.filtered <- filter(df, adj.pval < 0.01, log2FC > 0.5)
features <- list(df.filtered [,1])

DefaultAssay(sc)="RNA"
ctrl_genes = 50 #important

sc = AddModuleScore(object = sc, features = features, name = sample, ctrl = ctrl_genes)
FeaturePlot(sc, features = paste0(sample, '1'),pt.size = 1,order = T,label = FALSE) + scale_colour_gradient2(low = 'darkblue', mid = 'lightgrey', high = 'red', midpoint =0.5)
ggsave(filename = paste0(sample,"04_ForteScore_f.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"04_ForteScore_f.svg"), width=10 , height = 10)
VlnPlot(sc, features = paste0(sample, '1'),  group.by="Annotation_Level_1",pt.size = 0, sort = "decreasing") +NoLegend() + coord_flip() +geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)
ggsave(filename = paste(sample,"05_ForteScore_V.jpeg"), width=4, height = 5.5)
ggsave(filename = paste(sample,"05_ForteScore_V.svg"), width=4, height = 5.5)

#ttest_calculations for Scores
Naba.df<-FetchData(object = sc, vars=c("Annotation_Level_1", paste0(sample,"1")))
sink(paste0("",sample,"_06_TTest_Results.txt"))
print(paste0("t-test for Spp1 Mac WT IRI vs Ifn Mac in ", sample))
print(t.test(Naba.df[which(Naba.df$Annotation_Level_1 =="Spp1 Mac"),2], Naba.df[which(Naba.df$Annotation_Level_1=="IFN Mac"),2], alternative = "two.sided"))
sink()



