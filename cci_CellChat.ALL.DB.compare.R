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





## setting I 
## ---------------------------------- # 
OUTDIR_common <- "~/Spp1MacFibroisis/cellchat_results/"

DBs_col <- "all.db"
rds_path <- paste0(OUTDIR_common, DBs_col, "/")
rds  <- Sys.glob(file.path(rds_path, "*.rds"))
path_of_interest <- list( c("SEMA3", "SPP1", "FN1"))






## to compare CC1 vs CC2
## ---------------------------------- #
L_match <- c(rds[1], rds[3])
L_name <- rep(c("CXCL4_KO_IRI", "WT_IRI"),1)


for ( i in c(1)) {

        ## set up 
	## ---------------------------------- #

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

	print (DB)

	sub_cell_type_cc1 = levels(obj.L[["cc1"]]@idents)[c(6,9)]
        sub_cell_type_cc2 = levels(obj.L[["cc2"]]@idents)[c(6,9)]

       
	print (path_of_interest[[i]])
        path_of_interest_out <- path_of_interest[[i]]	



        ## find common / union pathways existed in each rds 
        ## ---------------------------------- #
	pathways.show <- CC2@netP$pathways[ CC2@netP$pathways %in% CC1@netP$pathways ]
        pathways.union <- union(CC1@netP$pathways, CC2@netP$pathways)



        ## 1st PDF 
	## ---------------------------------- #       
 
        pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.pdf"), width = 10, height = 10) 
        par(mfrow = c(1,2), xpd=TRUE)

	gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
	gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
	print(gg1 + gg2)
	


	## ---------------------------------- #
        ## compute the maximum number of cells per cell group and 
        ## the maximum number of interactions (or interaction weights) 
        ## across all datasets.
	## ---------------------------------- #

        weight.max <- getMaxWeight(obj.L, attribute = c("idents","count"))
        par(mfrow = c(1,2), xpd=TRUE)
        try( for ( stim in c(1:length(obj.L))) {
                netVisual_circle(obj.L[[stim]]@net$count,
                                                weight.scale = T,
                                                label.edge= F,
                                                edge.weight.max = weight.max[2],
                                                edge.width.max = 10,
                                                title.name = paste0("Number of interactions - ", names(obj.L)[stim]))
                })


        dev.off()




	## ---------------------------------- #
    	## Compare the total number of interactions and interaction strength
    	## ---------------------------------- #
	gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    	gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    	print(gg1 + gg2)
    
    
	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.heatmap.pdf"), width = 15, height = 7.5)
	gg1 <- netVisual_heatmap(cellchat)
	gg2 <- netVisual_heatmap(cellchat, measure = "weight")
	print(gg1 + gg2)
	dev.off()	
    	
    	




	## ---------------------------------- #
	## Compare the major sources and targets in 2D space
	## ---------------------------------- #
	## 2D : PDF 1 
	## ---------------------------------- #
	num.link <- sapply(obj.L, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.2D.pdf"), width = 8, height = 4)
	weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
	gg <- list()
	for ( d1 in 1:length(obj.L)) {
		  ##obj.L[[d1]] <- netAnalysis_computeCentrality(obj.L[[d1]], slot.name = "netP")
		  try(gg[[d1]] <- netAnalysis_signalingRole_scatter(obj.L[[d1]],
								    #label.size = 2, 
								    dot.size = c(1, 5),
								    title = names(obj.L)[d1], 
								    weight.MinMax = weight.MinMax))
	}
	try(print(patchwork::wrap_plots(plots = gg)))
	dev.off()





	## ---------------------------------- #
	## 2D : PDF 2 
	## ---------------------------------- #	 
	
	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.2D.selected.path.pdf"), width = 8, height = 4)
	gg <- list()
	


        for ( d2 in 1:length(obj.L)) {
                  ##obj.L[[d2]] <- netAnalysis_computeCentrality(obj.L[[d2]], slot.name = "netP")
		  try(gg[[d2]] <- netAnalysis_signalingRole_scatter(obj.L[[d2]],
								    #label.size = 2, 
								    title = names(obj.L)[d2],
								    dot.size = c(1, 5), 
								    weight.MinMax = weight.MinMax, 
								    signaling = path_of_interest_out))	
		  }
	try(print(patchwork::wrap_plots(plots = gg)))
	dev.off()

	


	## ---------------------------------- #
	## 2D : PDF 3,4 
	## ---------------------------------- #


	for ( sub_path in path_of_interest_out ) {
		pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.2D.selected.path.to.", sub_path, ".pdf"), width = 8, height = 4)
	        gg <- list()
	        for ( d3 in 1:length(obj.L)) {
	                  try(gg[[d3]] <- netAnalysis_signalingRole_scatter(obj.L[[d3]],
									    #label.size = 2, 
									    title = names(obj.L)[d3], 
									    dot.size = c(1, 5),
									    weight.MinMax = weight.MinMax, 
									    signaling = sub_path))
	                  }
	        try(print(patchwork::wrap_plots(plots = gg)))
	        dev.off()
		}
	





	## ---------------------------------- #
	## 2D : PDF 5 
	## ---------------------------------- #

	gg <- list()
	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.2D.MAC.FIBRO.pdf"), width = 11, height = 4)
	try(gg1 <- netAnalysis_signalingChanges_scatter(cellchat,  dot.size = 5, label.size=3, idents.use = "Mac"))
	try(gg2 <- netAnalysis_signalingChanges_scatter(cellchat,  dot.size = 5, label.size=3, idents.use = "Fibro"))
	try(print(patchwork::wrap_plots(plots = list(gg1,gg2))))
	dev.off()






        ## ---------------------------------- #
	## 2nd PDF 
	## ---------------------------------- #
	## To simplify the complicated network and gain insights
	## into the cell-cell communication at the cell type level,
	## we can aggregate the cell-cell communication
	## based on the defined cell groups.
	## ---------------------------------- #
	
	weight.max <- getMaxWeight(obj.L, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))

	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.netowrk..pdf"),, width = 10, height = 10) 
	
	par(mfrow = c(2,2), xpd=TRUE)
	for (sub_pathways.show in  pathways.show  ) { # pathways.union makes error 
		
        	try(weight.max <- getMaxWeight(obj.L, slot.name = c("netP"), attribute = sub_pathways.show)) # control the edge weights across different datasets
        	try(
        	for ( stim in 1:length(obj.L)) {
        	    netVisual_aggregate(obj.L[[stim]], signaling = sub_pathways.show, 
                                                                layout = "circle", 
                                                                edge.weight.max = weight.max[1], 
                                                                edge.width.max = 10, 
                                                                signaling.name = paste(sub_pathways.show, names(obj.L)[stim]))
        	    netVisual_aggregate(obj.L[[stim]], signaling = sub_pathways.show, 
                                                                layout = "chord", 
                                                                signaling.name = paste(sub_pathways.show, names(obj.L)[stim]))


		}
        	)}
    
    	dev.off()




	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.netowrk.II.pdf"), width = 10, height = 10)

	par(mfrow = c(2,2), xpd=TRUE)
        for (sub_pathways.show in  pathways.show  ) { # pathways.union makes error 

                try(weight.max <- getMaxWeight(obj.L, slot.name = c("netP"), attribute = sub_pathways.show)) # control the edge weights across different datasets
                try(num_ver_receiver <- obj.L[["cc1"]]@meta$labels %>% unique() %>% length())
                try(num_ver_receiver_half <- round(num_ver_receiver/2))
                vertex.receiver = seq(1,num_ver_receiver_half) #

                try(weight.max <- getMaxWeight(obj.L, slot.name = c("netP"), attribute = sub_pathways.show)) # control the edge weights across different datasets
                try(
                for ( stim in 1:length(obj.L)) {
                    netVisual_aggregate(obj.L[[stim]], signaling = sub_pathways.show,
                                                                 layout = "hierarchy",
                                                                 vertex.receiver = vertex.receiver,
                                                                 edge.weight.max = weight.max[1],
                                                                 edge.width.max = 10, 
                                                                 signaling.name = paste(sub_pathways.show, names(obj.L)[stim]))



                }
                )}

        dev.off()	



        ## ---------------------------------- #
	## 3rd PDF 
	## ---------------------------------- #


	for ( stim in 1:length(obj.L)) {
      	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.network", names(obj.L)[stim], ".pdf"), width = 10, height = 10) 
    
		for (sub_pathways.show in  pathways.show  ) { # pathways.union makes error 
		
			try(weight.max <- getMaxWeight(obj.L, slot.name = c("netP"), attribute = sub_pathways.show)) # control the edge weights across different datasets
			try(num_ver_receiver <- obj.L[[stim]]@meta$labels %>% unique() %>% length())
			try(num_ver_receiver_half <- round(num_ver_receiver/2))
			vertex.receiver = seq(1,num_ver_receiver_half) #
        			
			try(netVisual_aggregate(obj.L[[stim]], signaling = sub_pathways.show, 
                        	                          vertex.receiver = vertex.receiver, 
                        	                          edge.weight.max = weight.max[1], 
                        	                          edge.width.max = 10, signaling.name = paste(sub_pathways.show, names(obj.L)[stim])))
              		
			par(mfrow = c(3,2), xpd=TRUE)
			try(pairLR.path <- extractEnrichedLR(cellchat, signaling = sub_pathways.show, geneLR.return = FALSE))
			try(LR.show_1 <- pairLR.path[1:3,] %>% .[!is.na(.)] )
			try(LR.show_2 <- pairLR.path[4:6,] %>% .[!is.na(.)] )
			try(LR.show_3 <- pairLR.path[7:9,] %>% .[!is.na(.)] )
			try(LR.show_4 <- pairLR.path[10:12,] %>% .[!is.na(.)] )	
			try(LR.show_5 <- pairLR.path[13:15,] %>% .[!is.na(.)] )
			try(LR.show_6 <- pairLR.path[16:18,] %>% .[!is.na(.)] )
			try(LR.show_7 <- pairLR.path[19:21,] %>% .[!is.na(.)] )
			try(LR.show_8 <- pairLR.path[22:24,] %>% .[!is.na(.)] )
			try(LR.show_9 <- pairLR.path[25:27,] %>% .[!is.na(.)] )
			try(LR.show_10 <- pairLR.path[28:30,] %>% .[!is.na(.)] )
                    	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_1, vertex.receiver = vertex.receiver))
                    	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_2, vertex.receiver = vertex.receiver))
		      	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_3, vertex.receiver = vertex.receiver))	
			try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_4, vertex.receiver = vertex.receiver))
			try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_5, vertex.receiver = vertex.receiver))
			try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_6, vertex.receiver = vertex.receiver))
			try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_7, vertex.receiver = vertex.receiver))
			try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_8, vertex.receiver = vertex.receiver))
			try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_9, vertex.receiver = vertex.receiver))
			try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy", pairLR.use = LR.show_10, vertex.receiver = vertex.receiver))

			}
	dev.off()
	}	
  
	

	## ---------------------------------- #
        ## 4th PDF
	## ---------------------------------- #  
    
	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.network.bubble.pdf"), width = 10, height = 10) 
	par(mfrow = c(2,2), xpd=TRUE)


	#sub_cell_type_cc1 = levels(obj.L[["cc1"]]@idents)
	#sub_cell_type_cc1 = levels(obj.L[["cc2"]]@idents)

        try(gg1 <- netVisual_bubble(obj.L[["cc1"]], sources.use = sub_cell_type_cc1, 
                         remove.isolate = FALSE, 
                         angle.x = 90, 
                         font.size = 5,
                         title.name = paste0("as.sender.","cc1.lineage.cells")
                          ))
        try(gg1.cc2 <- netVisual_bubble(obj.L[["cc2"]], sources.use = sub_cell_type_cc2, 
                        remove.isolate = FALSE, 
                        angle.x = 90, 
                        font.size = 5,
                        title.name = paste0("as.sender.", "cc2.lineage.cells")
                          ))
        print (gg1+gg1.cc2)
        try(gg2 <- netVisual_bubble(obj.L[["cc1"]], targets.use = sub_cell_type_cc1, 
                         remove.isolate = FALSE, 
                         angle.x = 90,
                         font.size = 5,
                         title.name = paste0( ".as.targeted.", "cc1.lineage.cells")
                          ))
        try(gg2.cc2 <- netVisual_bubble(obj.L[["cc2"]], targets.use = sub_cell_type_cc2, 
                        remove.isolate = FALSE, 
                        angle.x = 90,
                        font.size = 5,
                        title.name = paste0( ".as.targeted.", "cc2.lineage.cells")
                          ))
        
        print (gg2+gg2.cc2)
        
       try(print(netVisual_bubble(obj.L[["cc1"]], 
                         remove.isolate = FALSE, 
                         angle.x = 90,
                         font.size = 5,
                         title.name = paste0( "cc1.all.cells"))
        ))
        try(print(netVisual_bubble(obj.L[["cc2"]], 
                         remove.isolate = FALSE, 
                         angle.x = 90,
                         font.size = 5,
                         title.name = paste0( "cc2.all.cells"))
        ))
        
        
        
        
        
    dev.off()



}


