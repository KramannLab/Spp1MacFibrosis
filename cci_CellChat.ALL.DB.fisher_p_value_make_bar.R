#!/usr/local/bin/Rscript
setwd("./")




## Written by Hyojin Kim
## Feb. 2022
## ---------------------------------- ##

library(Seurat)
library(dplyr)
library(cowplot)
library(tidyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(gtools)
library(ggforce)
library(furrr)
library(readr)
library(parallel)
library(tibble)
library(purrr)
library(reticulate)
library(phateR)
library(viridis)
library(ggrepel)
library(tidyr)
library(CellChat)
library(patchwork)
library(igraph)
library(cowplot)
library(stringr)
library(base)


args <- commandArgs(trailingOnly=TRUE)






## function 
## ------------------------------------------------- ## 
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









## setting I 
## mk dir and define filename
## ------------------------------------------------- ##
OUTDIR_common <- "~/Spp1MacFibroisis/cellchat_results/"

DBs_col <- "all.db"
rds_path <- paste0(OUTDIR_common, DBs_col, "/")
rds  <- Sys.glob(file.path(rds_path, "*.rds"))
path_of_interest <- list( c("SEMA3", "SPP1", "FN1"))










## To compare CC1 vs CC2
## ------------------------------------------------- ##
L_match <- c(rds[1], rds[3])
L_name <- rep(c("CXCL4_KO_IRI", "WT_IRI"),1)



for ( i in c(1)) {


	CC1 <- readRDS(paste0(L_match[i]))
	CC2 <- readRDS(paste0(L_match[i+1]))


        CC1 <- netAnalysis_computeCentrality(CC1, slot.name = "netP")
        CC2 <- netAnalysis_computeCentrality(CC2, slot.name = "netP")
	
	obj.L <- list(cc1 = CC1, cc2 = CC2)
	cellchat <- mergeCellChat(obj.L, add.names = names(obj.L), cell.prefix = TRUE)
	
	CC1_name <- L_name[i]
	CC2_name <- L_name[i+1]
	
	print (CC1_name)
	print (CC2_name)
	
	DB <- L_match[i] %>% strsplit(., ".cellcat.") %>% unlist() %>% .[2] %>% strsplit(., ".rds") %>% unlist() %>% .[1]
	rds_path_a <- L_match[i] %>% strsplit(., "[/]") %>% unlist() %>% .[1:7]
	rds_path <- paste(rds_path_a, collapse="/")
	
	
	make_bar <- function(cc1_mat, cc2_mat, key1, filename) {

		cc1_mat %>% filter(Var1== key1 | Var2== key1) -> cc1_key
		cc2_mat %>% filter(Var1== key1 | Var2== key1) -> cc2_key 
			
		cols <- c("Var1", "Var2")
		for ( col in cols ) {
			cc1_key[[col]] %>% as.character() -> cc1_key[[col]]
			cc2_key[[col]] %>% as.character() -> cc2_key[[col]]
			}
		
		
		cc1_key %>% mutate(key = paste0(pmin(Var1, Var2), "-", pmax(Var1, Var2))) %>% 
			group_by(key) %>% 
			summarise(sum = sum(value)) %>% 
			arrange(desc(sum)) %>% 
	       		mutate(perc = sum / sum(sum)) -> cc1_key_count
		
		cc1_key_count$group <- "cc1"
		
		cc2_key %>% mutate(key = paste0(pmin(Var1, Var2), "-", pmax(Var1, Var2))) %>% 
	       	        group_by(key) %>% 
	                summarise(sum = sum(value)) %>%
	                arrange(desc(sum)) %>% 
			mutate(perc = sum / sum(sum)) -> cc2_key_count
		
		cc2_key_count$group <- "cc2"
		
		cc_key <- rbind(cc1_key_count, cc2_key_count)
		word_rp <- paste0(key1, "-|-", key1)
		cc_key$key %>% gsub(word_rp, "", .) -> cc_key$new_key 
		cc_key$new_key %>% factor(., levels=sort(unique(cc_key$new_key), decreasing = TRUE)) -> cc_key$new_key
		
		pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".", key1, ".", filename, ".barplot.pdf"), width = 6, height = 6) 
		p <- ggplot(data = cc_key, aes(x = sum, y = new_key)) + 
			  geom_bar(stat = "identity", 
		           aes(fill = factor(group)), 
		           position = "dodge") +
			   theme_minimal() +
			   labs(title="", 
			        x="", 
			        y= paste0("The ", filename, " of interactions with ", key1))
		print(p)
		dev.off()

		return(cc_key)
		}

	
	obj.L$cc1@net$count %>% melt() -> cc1_mat
	obj.L$cc2@net$count %>% melt() -> cc2_mat
	key1 <- "Mac"
	cc_mat_out <- make_bar(cc1_mat, cc2_mat, key1, "count")


	obj.L$cc1@net$weight %>% melt() -> cc1_mat_w
	obj.L$cc2@net$weight %>% melt() -> cc2_mat_w
	make_bar(cc1_mat_w, cc2_mat_w, key1, "weight")





	## run fisher exact test 
        ## ------------------------------------------------- ##
	cc_mat_out1 <- cc_mat_out %>% filter(group=="cc1") %>% select(c(-perc, -key))
	cc_mat_out2 <- cc_mat_out %>% filter(group=="cc2") %>% select(c(-perc, -key))
	cc_fisher <- inner_join(cc_mat_out1, cc_mat_out2, by="new_key") %>% select(new_key, sum.x, sum.y) %>% as.data.frame()
	colnames(cc_fisher) <- c("receiver", CC1_name, CC2_name)
	cc_fisher <- cc_fisher %>% column_to_rownames("receiver")
	
	OUTDIR <- rds_path
        cc_fisher <- as.matrix(cc_fisher)	
	cc_fisher_out <- run_fisher(cc_fisher, "cellchat.count.")





	## differential  R L 
	## ------------------------------------------------- ##	
	
        sub_cell_type_cc1 = levels(obj.L[["cc1"]]@idents)[c(6,9)]
        sub_cell_type_cc2 = levels(obj.L[["cc2"]]@idents)[c(6,9)]

	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".diff.heatmap.MAC.FIBRO.SEMA3.pdf"), width = 4, height = 4)
	gg1 <- netVisual_bubble(cellchat, sources.use = sub_cell_type_cc1, 
						signaling = path_of_interest[[1]],
						targets.use = sub_cell_type_cc1,  
						comparison = c(1, 2),
						angle.x = 90, remove.isolate = T)
	print(gg1)
	dev.off()


        pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".diff.heatmap.MAC.FIBRO.allpath.pdf"), width = 5, height = 12)
        gg1 <- netVisual_bubble(cellchat, sources.use = sub_cell_type_cc1,
                                                ##signaling = path_of_interest[1],
                                                targets.use = sub_cell_type_cc1,
                                                comparison = c(1, 2),
                                                angle.x = 90, remove.isolate = T)
        print(gg1)
        dev.off()




	}# close main loop
