#!/usr/local/bin/Rscript
setwd("./")


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


args <- commandArgs(trailingOnly=TRUE)





## Written by Hyojin Kim 
## Feb 2022
## ---------------------------------- #





## setting I 
## mk dir and define filename
## ---------------------------------- # 

OUTDIR_common <- "~/Spp1MacFibroisis/cellchat_results/"
OUTDIR <- c("~/Spp1MacFibroisis/cellchat_results/all.db/") 

if(!dir.exists(OUTDIR_common)) dir.create(OUTDIR_common)
if(!dir.exists(OUTDIR)) dir.create(OUTDIR)


rds_path <- ""~/Spp1MacFibroisis/"
rds <- c("WTvsCxcl4KO_IRI_snRNA_Integrated.h5ad.to.RDS") 
DBs_col <- "all.db"


rds1 <- "leukocytes_cxcl4_subclusters_WTvsCxcl4KO_IRI_snRNA_Integrated.h5ad.to.RDS" 









## setting II
## add "cell_type" column same as the cell type meta of your object 
# ---------------------------------- #

remove_celltype = c("Leukocytes", "Adipocytes") # if nothing to remove, just type "nothing"
add_celltype = c("Mac1", "Mac2", "Bcells", "Tcells")











## To add "cell_type" column same as the cell type meta of your object
## To make subsets for each condition 
## ---------------------------------- #

R <- readRDS(paste0( rds_path, rds ))
# to remove ":_" in cell type name 
R$cell_type <- R$leiden0.5 %>% as.character() %>% strsplit(., "[: ]") %>% lapply(., function(x) x[3]) %>% unlist()
assay_name <- DefaultAssay(R)
conditions <- R$condition %>% unique() %>% as.character() %>% mixedsort()
R_orig <- R


L <- readRDS(paste0( rds_path, rds1 ))
L <- subset(L, leiden0.2 %in% add_celltype)
L$leiden0.2 %>% as.character() -> L$leiden0.2
Idents(L) <- "leiden0.2"









## merge
## ---------------------------------- #

R <- R[,!R$cell_type %in% remove_celltype ]
RL <- merge(R, y = L, add.cell.ids = c("all", "leuko"), project = "cxcl4")
RL$cell_type_origin <- RL$cell_type
RL$leiden0.2 %>% as.character() -> RL$leiden0.2
RL$leiden0.2 %>% unique() %>% .[2:5] -> valid_list
RL$cell_type <- RL[[]] %>% mutate(cell_type = ifelse(leiden0.2 %in% valid_list, leiden0.2, cell_type_origin)) %>% select(cell_type) 
RL$condition %>% toupper() -> RL$condition
RL$condition %>% unique() %>% sort() -> conditions
RL$cell_type %>% gsub("Mac2|Mac1", "Mac", .) -> RL$cell_type


print ('merge done')









## To simplify conditions'name 
## ---------------------------------- #
## > conditions
## [1] "CXCL4_KO_IRI"  "Cxcl4_KO_Sham" "WT_IRI"        "WT_Sham" 
## > conditions_SS
## [1] "C.I." "C.S." "W.I." "W.S."
## ---------------------------------- #

conditions %>% str_extract_all(., "\\b[A-Za-z]+") %>% toupper() %>% str_extract_all(., "\\b[A-Za-z]") -> conditions_S1
conditions %>% str_extract_all(., "[A-Za-z]+\\b") %>% toupper() %>% str_extract_all(., "\\b[A-Za-z]") -> conditions_S2
conditions_SS <- paste0(conditions_S1, ".", conditions_S2)

rds_col <- c()
for ( condition_a in c(conditions)) {
	rds_col[[condition_a]] <- subset(RL, condition == condition_a )
	}



for ( i in c(1:length(rds_col))) {

	S <- rds_col[[i]]
	mat <- S[[assay_name]]@data %>% as.data.frame(.) %>% rownames_to_column(.)
	colnames(mat)[1] <- "Gene" 

	# to add prefix to the cell type or cell name 
	# meta_dat <- cbind( rownames(S@meta.data), paste0(conditions_SS[[i]], ".", S$cell_type )) %>% as.data.frame()
  	meta_dat <- cbind( rownames(S@meta.data), S$cell_type ) %>% as.data.frame()

	colnames(meta_dat) <- c("Cell","cell_type")
  
	write.table(meta_dat, file=paste0(OUTDIR_common, Project(S), ".", conditions[[i]], ".filtered.meta.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote = FALSE)
  
	  # cellchat starts 
	  # Part I: Data input & processing and initialization of CellChat object
	  # ---------------------------------- #
	  data.input = S[[assay_name]]@data %>% as.matrix() # normalized data matrix
	  colnames(meta_dat) <- c("Cell","labels")
	  meta <- meta_dat %>% as.data.frame()
	  meta$labels <- as.factor(meta$labels) %>% as.character()
	  rownames(meta) <- meta$Cell
	  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
	  cellchat <- addMeta(cellchat, meta = meta)
	  cellchat <- setIdent(cellchat, ident.use = "labels")
	  groupSize <- as.numeric(table(cellchat@idents))
  


		# Set the ligand-receptor interaction database
		# ---------------------------------- #
        	CellChatDB <- CellChatDB.mouse
	  	CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor") ) # [1] annotation =="Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
        	cellchat@DB <- CellChatDB.use


       		# Preprocessing the expression data 
		# for cell-cell communication analysis
        	# ---------------------------------- #
		cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
        	future::plan("multiprocess", workers = 10) # do parallel
        	cellchat <- identifyOverExpressedGenes(cellchat)
        	cellchat <- identifyOverExpressedInteractions(cellchat)
        	cellchat <- projectData(cellchat, PPI.mouse)
        
       		# Part II: Inference of cell-cell communication network
        	# Compute the communication probability and infer cellular communication network
		# ---------------------------------- #
        	try(cellchat <- computeCommunProb(cellchat, raw.use = TRUE))
        	try(cellchat <- filterCommunication(cellchat, min.cells = 10))

       		# Extract the inferred cellular communication network as a data frame
		# ---------------------------------- #
        	try(df.net <- subsetCommunication(cellchat))
        	
		# Infer the cell-cell communication at a signaling pathway level
        	# ---------------------------------- #
		try(cellchat <- computeCommunProbPathway(cellchat))

        	# Calculate the aggregated cell-cell communication network
		# ---------------------------------- #
       		try(cellchat <- aggregateNet(cellchat))
        	try(groupSize <- as.numeric(table(cellchat@idents)))
        
      
	       	pdf(paste0(OUTDIR, Project(S), ".", conditions[[i]],  ".cell.chat.agg.cell.cell.network.pdf"), width = 10, height = 10)
        
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
                          } # close for (n_mat in 1:nrow(mat))
                          )
                    
              
                # Part III: Visualization of cell-cell communication network
                # ---------------------------------- #
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
                        ggsave(filename=paste0(OUTDIR, Project(S), ".", conditions[[i]], ".",  n, "_L-R_contribution.12.Mar.pdf"), plot=gg, width = 3, height = 2, units = "in", dpi = 300)
                                        
                        }) # close for ( n in c(cellchat@netP$pathways) ){
            
        dev.off()
        
        try(saveRDS(cellchat, paste0(OUTDIR, Project(S), ".", conditions[[i]],  ".cellcat.", DBs_col, ".rds")))
 
  

	cell_types <- cellchat@meta$labels %>% unique() %>% as.character() 
        cell_types_chr <- length(cell_types)

	pdf(paste0(OUTDIR, Project(S), ".", conditions[[i]], ".network.pdf"), width = cell_types_chr/2.5, height = 5.5)
	       for ( sub_cell_type in cell_types ){
                  try(print(netVisual_bubble(cellchat, sources.use = sub_cell_type, remove.isolate = FALSE)))
                  try(print(netVisual_bubble(cellchat, targets.use = sub_cell_type, remove.isolate = FALSE)))
                  }
        dev.off()



	

  }



