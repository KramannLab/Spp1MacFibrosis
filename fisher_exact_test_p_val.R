#!/usr/local/bin/Rscript
setwd("./")




## Written by Hyojin Kim 
## Feb. 2022 
## ---------------------------------- ##


library(Seurat)
library(dplyr)
library(reshape2)
library(tidyverse)
library(gtools)
library(furrr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(stats)
library(rstatix)

args <- commandArgs(trailingOnly=TRUE)






## Setting
## ---------------------------------- ##

OUTDIR <- c("~/Spp1MacFibroisis/fisher_exact_test/") 
if(!dir.exists(OUTDIR)) dir.create(OUTDIR)







## Read file 
## ---------------------------------- ##
## "WTvsCxcl4KO_IRI_snRNA_Integrated.h5ad" -> "WTvsCxcl4KO_IRI_snRNA_Integrated.h5ad.to.RDS"
## refer to : https://github.com/genehaus/singlecellRNA-seq_2022/blob/main/scanpy_to_seurat/h5ad_to_rds.R
## ---------------------------------- ##

rds_path <- "~/Spp1MacFibroisis/fisher_exact_test/"
irds <- args[1] # "WTvsCxcl4KO_IRI_snRNA_Integrated.h5ad.to.RDS"
iremove_celltype = c("Adipocytes")








## To add "cell_type" column same as the cell type meta of your object
## To make subsets for each condition 
## ---------------------------------- ##

R <- readRDS(paste0( rds_path, rds ))

if ( rds == "WTvsCxcl4KO_IRI_snRNA_Integrated.h5ad.to.RDS")  {

		R$cell_type <- R$leiden0.5 %>% as.character() %>% strsplit(., "[: ]") %>% lapply(., function(x) x[3]) %>% unlist()
		assay_name <- DefaultAssay(R)
		conditions <- R$condition %>% unique() %>% as.character() %>% mixedsort()
		R <- R[,!R$cell_type %in% remove_celltype ]

	} else if ( rds == "Rao_MPC_07_seurat_obj_after_Clusterexclusion.rds" ) { ## public data 

		assay_name <- DefaultAssay(R)
		R$cell_type <- R$seurat_clusters %>% as.character()
	        R$condition <- R$disease	
		conditions <- R$disease %>% unique() %>% as.character() %>% mixedsort()

	} else if ( rds == "huCKD_CD10Neg_MPC_05_seurat_object.RDS" ) { ## public data
 
		Idents(R) <- "Annotation_Level_1"
		R$condition <- R$Kidney.Function
		R$cell_type <- Idents(R) %>% as.character()
		assay_name <- DefaultAssay(R)
	}










## Perform fisher-exact-test 
## ---------------------------------- #
## Refer to : https://www.reneshbedre.com/blog/fisher-exact-test.html
## ---------------------------------- #

ALL <- table(R$cell_type, R$condition) 
ALL_df <- ALL %>% as.data.frame() %>% 
			dcast(., Var1~Var2, value.var="Freq") %>% 
			column_to_rownames(var="Var1") %>% 
			as.matrix()

if ( rds == "scvi_output_clustered.h5ad.to.RDS")  {
		ALL_df[,c(4,3)] -> WT_df
		ALL_df[,c(2,1)] -> KO_df
		ALL_df[,c(2,4)] -> KO_WT_IRI
	} else if ( rds == "Rao_MPC_07_seurat_obj_after_Clusterexclusion.rds"  ) { 
		#ALL_df[,c(1,2)] -> DCM
		#ALL_df[,c(3,2)] -> ICM
		ALL_df <- as.data.frame(ALL_df)
		ALL_df$DICM <- as.numeric(ALL_df$DCM) + as.numeric(ALL_df$ICM)
		ALL_df[,c(4,2)] -> DICM
	} else if ( rds == "huCKD_CD10Neg_MPC_05_seurat_object.RDS" ) {
		ALL_df -> public_df 
	}


run_fisher <- function(raw, key) {
		out <- c()
		for ( n in seq(1,length(rownames(raw)))) {
			subject <- raw[n,]
			others <- colSums(raw) - subject 
			df <- rbind(subject, others)
			rownames(df)[1] <- rownames(raw)[n]
			out[[n]] <- fisher_test(as.matrix(df))
			out[[n]]$name <- rownames(raw)[n]
			}
		out_df <- out %>% do.call("rbind", .) %>% data.frame()
		out_df$adj.p <- p.adjust(out_df$p, method = "fdr", n = length(out_df$p))
		out_df$adj.p.significance <- stars.pval(out_df$adj.p)
		write.table(out_df, file=paste0(OUTDIR, "fisher_exact_test.fdr.", key, ".txt"), row.names=TRUE, col.names=TRUE, sep="\t")
		return (out_df)
		}



if ( rds == "scvi_output_clustered.h5ad.to.RDS")  {
		WT_df_out <- run_fisher(WT_df, "WT_df")
		KO_df_out <- run_fisher(KO_df, "KO_df")
		KO_WT_IRI_out <- run_fisher(KO_WT_IRI, "KO_WT_IRI")
	} else if (rds == "Rao_MPC_07_seurat_obj_after_Clusterexclusion.rds" ) {
		DICM_out = run_fisher(DICM, "DICM_df")
	} else if ( rds == "huCKD_CD10Neg_MPC_05_seurat_object.RDS" ) {
		public_df_out <- run_fisher(public_df, "public_df")
		}	





