# single cell RNA seq data analysis <br> to study intercellular communications with cellchat <br> Written by H.J. Kim <br>  In Jan 2022

<br>
<br>


from Seurat object to cellchat <br>  
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
	
	if(!dir.exists(OUTDIR_common)) dir.create(OUTDIR_common)
	OUTDIR %>% lapply(., function(x) if(!dir.exists(x)) dir.create(x))





#### 2-2. Read RDS data 



	rds_path <- "/xxx/"
	rds <- c("Seurat_object.RDS")
	R <- readRDS(paste0( rds_path, rds ))




#### 2-3. Make some list



	DBs_col <- c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")
	remove_celltype = c("Leukocytes", "Adipocytes") # if nothing to remove, just type "nothing" : remove_celltype = c("nothing")



#### 2-4. 
#### Make "cell_type" column of meta of your Seurat object which saves the cell annotation 
#### and "condition" for saving batch or disease of interest you want to divide when running cellchat



	R$cell_type <- R$column_name_save_cell_annotation 
	R$condition <- R$column_name_save_group_information
	
	# R[[]] %>% head() 
	# - + - + - + - + - + - + - + - + - + - + - + - #
	# rownames  # cell_type # condition # ... # ... # 
	# - + - + - + - + - + - + - + - + - + - + - + - #
	# AAACCCA.. # T cell    # disease   # ... # ... #
	# AAACCCG.. # T cell    # disease   # ... # ... #
	# AAACCCT.. # B cell    # disease   # ... # ... #
	# AAACCCC.. # B cell    # normal    # ... # ... #
	# ........  # T cell    # normal    # ... # ... #


	# if your data is coming from scanpy 
	R$cell_type <- R$cell_type %>% as.character() %>% strsplit(., "[: ]") %>% lapply(., function(x) x[3]) %>% unlist()
	
	# In the case of the DefaultAssay, 
	# "RNA" from Seurat object is ideal ("originalexp" from Scanpy)
	# Please do not use the one from integrated assay
	assay_name <- DefaultAssay(R)

	print (assay_name) # to check if the assay is "RNA" or "originalexp"
	 
	# 'condition' column has the exp. condition information in this R object
	conditions <- R$condition %>% unique() %>% as.character() %>% mixedsort()




#### 2-5. Simplify conditions' name


	# option 1 
	conditions %>% str_extract_all(., "\\b[A-Za-z]+") %>% toupper() %>% str_extract_all(., "\\b[A-Za-z]") -> conditions_S1
	conditions %>% str_extract_all(., "[A-Za-z]+\\b") %>% toupper() %>% str_extract_all(., "\\b[A-Za-z]") -> conditions_S2
	conditions_SS <- paste0(conditions_S1, ".", conditions_S2)
	
	#> conditions
	#[1] "ABC_CELL_HEALTHY"  "ABC_CELL_DISEASE" "DEF_CELL_HEALTHY" "DEF_CELL_DISEASE"
	#> conditions_SS
	#[1] "A.H." "A.D." "D.H." "D.D."
	
	# option 2 (typing by manual)
	# it should be the same number of the conditions in 2-4.
	conditions_SS <- c("A.H.", "A.D.", "D.H.", "D.D.")



#### 2-6. Subset R object based on the condition 




	rds_col <- c()
	for ( condition_a in c(conditions)) {
        	R <- R[,!R$cell_type %in% remove_celltype ]
        	rds_col[[condition_a]] <- subset(R, condition == condition_a )
	}




#### 2-7. Run CellChat 


The goal is to generate pdf files which inclue all possible cellchat output for each cell type and each pathway in a given condition and database of cellchat



	for ( i in c(1:length(rds_col))) {
	
	        S <- rds_col[[i]]
	        mat <- S[[assay_name]]@data %>% as.data.frame(.) %>% rownames_to_column(.)
	        colnames(mat)[1] <- "Gene"
	
	        # to add prefix to the cell type with the output from step 2-5. 
	        meta_dat <- cbind( rownames(S@meta.data), paste0(conditions_SS[[i]], ".", S$cell_type )) %>% as.data.frame()
	
	        colnames(meta_dat) <- c("Cell","cell_type")
	        write.table(meta_dat, file=paste0(OUTDIR_common, Project(S), ".", conditions[[i]], ".filtered.meta.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote = FALSE)
		
		
		# cellchat starts 
		
	          # Part I: Data input & processing and initialization of CellChat object
		
		
	          data.input = S[[assay_name]]@data %>% as.matrix() # normalized data matrix
	          colnames(meta_dat) <- c("Cell","labels")
	          meta <- meta_dat %>% as.data.frame()
	          meta$labels <- as.factor(meta$labels) %>% as.character()
	          rownames(meta) <- meta$Cell
	          cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
        	  cellchat <- addMeta(cellchat, meta = meta)
        	  cellchat <- setIdent(cellchat, ident.use = "labels")
        	  groupSize <- as.numeric(table(cellchat@idents))
		
		
		  for ( sub_db in c(1:length(DBs_col)) ) {
	
	                # Set the ligand-receptor interaction database
	                CellChatDB <- CellChatDB.mouse # if human, type CellChatDB.human
	                CellChatDB.use <- subsetDB(CellChatDB, search = DBs_col[sub_db] ) # [1] annotation =="Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
		        cellchat@DB <- CellChatDB.use
	
	                # Preprocessing the expression data for cell-cell communication analysis
	                cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
	                future::plan("multiprocess", workers = 4) # do parallel
		                cellchat <- identifyOverExpressedGenes(cellchat)
        	        cellchat <- identifyOverExpressedInteractions(cellchat)
        	        cellchat <- projectData(cellchat, PPI.mouse)
	
		
	
	                # Part II: Inference of cell-cell communication network
	
		
        	        # Compute the communication probability and infer cellular communication network
        	        try(cellchat <- computeCommunProb(cellchat, raw.use = TRUE))
        	        try(cellchat <- filterCommunication(cellchat, min.cells = 10))
        	        # Extract the inferred cellular communication network as a data frame
                	try(df.net <- subsetCommunication(cellchat))
                	# Infer the cell-cell communication at a signaling pathway level
         	       try(cellchat <- computeCommunProbPathway(cellchat))
                	# Calculate the aggregated cell-cell communication network
                	try(cellchat <- aggregateNet(cellchat))
                	try(groupSize <- as.numeric(table(cellchat@idents)))

	
		pdf(paste0(OUTDIR[sub_db], Project(S), ".", conditions[[i]],  ".cell.chat.agg.cell.cell.network.", DBs_col[sub_db], ".pdf"), width = 10, height = 10)
		
        	        par(mfrow = c(1,2), xpd=TRUE)
        	        try(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions"))
        	        try(netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength"))
	
		                try(mat <- cellchat@net$weight)
        	        try(par(mfrow = c(3,3), xpd=TRUE))
        		        #par(mar=c(1,1,1,1))
                	try(
                	for (n_mat in 1:nrow(mat)) {
                	          mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
                	          mat2[n_mat, ] <- mat[n_mat, ]
                	          netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[n_mat])
                	          }
                	          )
		
		
                	# Part III: Visualization of cell-cell communication network
		
				
                	try(num_ver_receiver <- meta$labels %>% unique() %>% length())
                	try(num_ver_receiver_half <- round(num_ver_receiver/2))
                	try(
                	for ( n in c(cellchat@netP$pathways) ){
	
                        # ANGPTL, ncWNT, MK, PTN
                        #par(mfrow = c(3,2), xpd=TRUE)
                        pathways.show <- n
                        vertex.receiver = seq(1,num_ver_receiver_half) # a numeric vector. 
                        try(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy", vertex.receiver = vertex.receiver))
                        try(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle"))
                        try(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord"))
                        try(netAnalysis_contribution(cellchat, signaling = pathways.show))
                        try(pairLR.path <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE))
                        try(path_legnth <- dim(pairLR.path)[1])
                        try(print (path_legnth))
                        try(LR.show <- pairLR.path[,])
                        for ( sub_LR in LR.show ) {
                            try(netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = sub_LR, layout = "hierarchy", vertex.receiver = vertex.receiver))
                            try(netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = sub_LR, layout = "circle"))
                            try(netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = sub_LR, layout = "chord"))

                            }

                        try(gg <- netAnalysis_contribution(cellchat, signaling = pathways.show ))
                        ggsave(filename=paste0(OUTDIR[sub_db], Project(S), ".", conditions[[i]], ".",  n, "_L-R_contribution.12.Mar.pdf"), plot=gg, width = 3, height = 2, units = "in", dpi = 300)

			})
	
		dev.off() # close the pdf 
		
		
		try(saveRDS(cellchat, paste0(OUTDIR[sub_db], Project(S), ".", conditions[[i]],  ".cellcat.", DBs_col[sub_db], ".rds")))
		
		
		
        	cell_types <- cellchat@meta$labels %>% unique() %>% as.character()
        	cell_types_chr <- length(cell_types)
		
	        pdf(paste0(OUTDIR[sub_db], Project(S), ".", conditions[[i]], ".network.pdf"), width = cell_types_chr/2.5, height = 5.5)
        	       for ( sub_cell_type in cell_types ){
        	          try(print(netVisual_bubble(cellchat, sources.use = sub_cell_type, remove.isolate = FALSE)))
        	          try(print(netVisual_bubble(cellchat, targets.use = sub_cell_type, remove.isolate = FALSE)))
        	          }
        	dev.off()
	
	
		}}

	 

<br>
<br>



### 3. acknowledgement

	
#### 3-1. website
	
	https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html



<br>
<br>




### 6. SessionInfo



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






