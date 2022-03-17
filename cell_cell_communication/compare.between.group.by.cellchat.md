# single cell RNA seq data analysis  <br> to study intercellular communications II with cellchat <br> (for comparison between groups) <br> Written by H.J. Kim <br>  In Jan 2022

<br>
<br>


Next step to README.md <br> 
To compare cellchat outputs between groups <br>  
This analysis is based on the data from mice.

### 1. Call libraries <br> 



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
	library(phateR)
	library(ggrepel)
	library(CellChat)
	library(patchwork)
	library(igraph)
	library(stringr)

	
	

### 2. Analyse scRNA-seq data from a single sample/condition/batch <br> 
Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html


### 2-1. Make output directories



	setwd("./") 
	
	OUTDIR_common <- "./cellchat_results/"
	OUTDIR <- c("./cellchat_results/Secreted_Signaling/",
	            "./cellchat_results/ECM-Receptor/",
	            "./cellchat_results/Cell-Cell_Contact/" )
	



#### 2-2. Read RDS data saved by README.md 

	
	DBs_col <- c("Secreted_Signaling","ECM-Receptor","Cell-Cell_Contact")
	rds_path <- paste0(OUTDIR_common, DBs_col, "/")
	rds  <- Sys.glob(file.path(rds_path, "*.rds"))
	
	
	# + - + - + - + - + - + - + - + - + - + - + - + - #
	# select the rds file to compare 
	# In example below, 
	# rds[1] vs rds[3]  
	# rds[5] vs rds[7]  
	# rds[9] vs rds[11] 
	# 
	# make "L_name" which represents the short name of conditions
	# rds[1], rds[5], rds[9] is from the condition "A.H"
	# rds[3], rds[7], rds[11] is from the condition "D.H" 
	# + - + - + - + - + - + - + - + - + - + - + - + - #
	
	L_match <- c(rds[1], rds[3], rds[5], rds[7], rds[9], rds[11])
	L_name <- rep(c("A.H", "D.H"),3)



#### 2-3. Run cellchat 


	# for each "loop" 
	# I plan to call 2 rds together to compare between those 2 rds in the list of "L_match" 
	# So, I set up the call_numbers with odd numbers  
	call_numbers <- c(1,3,5) 


	for ( i in call_numbers) {

        	# + - + - + - + - + - + - + - + - + - + - + - + - # 
        	# set up 
        	# + - + - + - + - + - + - + - + - + - + - + - + - #
		
		
        	CC1 <- readRDS(paste0(L_match[i]))
        	CC2 <- readRDS(paste0(L_match[i+1]))
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
		
        	sub_cell_type_cc1 = levels(obj.L[["cc1"]]@idents)[c(6,9,10)]
        	sub_cell_type_cc2 = levels(obj.L[["cc2"]]@idents)[c(6,9,10)]
		
		
		
        	# + - + - + - + - + - + - + - + - + - + - + - + - # 
        	# find common / union pathways existed in each rds 
        	# + - + - + - + - + - + - + - + - + - + - + - + - # 
        	pathways.show <- CC2@netP$pathways[ CC2@netP$pathways %in% CC1@netP$pathways ]
        	pathways.union <- union(CC1@netP$pathways, CC2@netP$pathways)
			
		
		
		
		
		# + - + - + - + - + - + - + - + - + - + - + - + - # 
        	# 1st PDF
		#
		# ref : https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
		# "Compare the total number of interactions and interaction strength"  
        	# + - + - + - + - + - + - + - + - + - + - + - + - #
		
        	pdf(paste0(rds_path, CC1_name, ".vs.", CC2_name, ".",  DB, ".compare.pdf"), width = 10, height = 10)
        	par(mfrow = c(1,2), xpd=TRUE)
			
		
		
        	# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + #  
        	# Compare the total number of interactions and interaction strength
        	# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + #  
        	gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
        	gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
        	print(gg1 + gg2)		
		
		
	        # + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + #
	        # compute the maximum number of cells per cell group and 
	        # the maximum number of interactions (or interaction weights) 
	        # across all datasets.
	        # + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + #
	
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
		
		 dev.off() # close 1st pdf 
		
		
		
		
		
		# + - + - + - + - + - + - + - + - + - + - + - + - # 
        	# 2nd PDF 
		#
		# ref : https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
		# "Visualize the inferred signaling network using the lifted object" at a pathway level 
        	# + - + - + - + - + - + - + - + - + - + - + - + - #
		
        	# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + #
        	# To simplify the complicated network and gain insights
        	# into the cell-cell communication at the cell type level,
        	# we can aggregate the cell-cell communication
        	# based on the defined cell groups.
        	# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + #
		
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
		
        	dev.off() # close the 2nd pdf 
		
		
		
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
		
        	dev.off() # close the II version of 2nd pdf 
		
		
		
		
		
        	# + - + - + - + - + - + - + - + - + - + - + - + - # 
	        # 3rd PDF 
		# 
		# ref : https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html	
		# "Visualize the inferred signaling network using the lifted object" at the single ligand and receptor level 
	        # + - + - + - + - + - + - + - + - + - + - + - + - #
		
		
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
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_1, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_2, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_3, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_4, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_5, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_6, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_7, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_8, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_9, vertex.receiver = vertex.receiver))
                        	try( netVisual_individual(obj.L[[stim]], signaling = sub_pathways.show, layout = "hierarchy", pairLR.use = LR.show_10, vertex.receiver = vertex.receiver))

                        	}
	        	 dev.off() # close 3rd pdf
       		 } # close for the close "for loop"
		
		
		
		# + - + - + - + - + - + - + - + - + - + - + - + - #
        	# 4th PDF
		# 
		# ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
		# "Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs"
        	# + - + - + - + - + - + - + - + - + - + - + - + - #
		
		
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
		
		
		dev.off() # close 4th pdf 

	} # close the top-level for loop



<br>



### 3. acknowledgement

	
#### 3-1. website
	
	https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
	https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html

<br>
<br>




### SessionInfo



```
## > sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 8

locale:
 [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C             
 [3] LC_TIME=en_US.utf8        LC_COLLATE=en_US.utf8    
 [5] LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8   
 [7] LC_PAPER=en_US.utf8       LC_NAME=C                
 [9] LC_ADDRESS=C              LC_TELEPHONE=C           
[11] LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C      

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] patchwork_1.1.1     CellChat_1.1.3      Biobase_2.52.0     
 [4] BiocGenerics_0.38.0 igraph_1.2.11       ggrepel_0.9.1      
 [7] viridis_0.6.2       viridisLite_0.4.0   phateR_1.0.7       
[10] Matrix_1.4-0        furrr_0.2.3         future_1.23.0      
[13] ggforce_0.3.3       gtools_3.9.2        forcats_0.5.1      
[16] stringr_1.4.0       purrr_0.3.4         readr_2.1.1        
[19] tibble_3.1.6        tidyverse_1.3.1     ggplot2_3.3.5      
[22] reshape2_1.4.4      tidyr_1.1.4         cowplot_1.1.1      
[25] dplyr_1.0.7         SeuratObject_4.0.4  Seurat_4.0.6
```






