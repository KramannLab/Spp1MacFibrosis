#!/usr/local/bin/Rscript

setwd("./")


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
library(monocle) # version 2.20
library(slingshot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(psych)
library(stats)
library(cartography)
library(colorRamps)
library(scales)
library(colorspace)
library(rcartocolor)
args <- commandArgs(trailingOnly=TRUE)
set.seed(42)
source("./monocle.function.R")








## Written by Hyojin Kim
## Feb. 2022
## ---------------------------------- ##







## Setting I 
## ---------------------------------- # 
outdir <- c("~/Spp1MacFibroisis/slingshot_out/")
if(!dir.exists(outdir)) dir.create(outdir)

rds_path <- "~/Spp1MacFibroisis/"
rds <- "MI_ImmuneCell_13_SeuratObjt_AfterFirstAnalysis.RDS" ## Publicly available data
cell_type <- "Annotation_Level_1"











## Read
## ---------------------------------- #
R <- readRDS(paste0( rds_path, rds ))

include_celltype <- c("Ly6c2hi Mono", "IFN Mac", "Arg1 Mac")
R$Annotation_Level_1 %>% gsub("/", ".", .) %>% as.character() -> R$Annotation_Level_1
Idents(R) <- "Annotation_Level_1"
R$cell_type <- Idents(R) %>% as.character()
R@assays$ECMRegulators=NULL
R_ALL <- R










## Colors
## ---------------------------------- #
real_colors <- c("Res Mac" = "#ee756e", 
	    	"Granulocytes" = "#cd9406", 
	    	"Arg1 Mac" = "#7baf2a", 
	    	"B Cells" = "#3aae66", 
	    	"Ly6c2hi Mono" = "#27b5bf", 
	    	"cDC2" = "#4498d4",
	    	"T-Cells.NK Cells" = "#a080b9", 
	    	"IFN Mac" = "#d76da8") 



umap_colors <- c("Res Mac" = "#D3D3D3",
                "Granulocytes" = "#D3D3D3",
                "Arg1 Mac" = "#7baf2a",
                "B Cells" = "#D3D3D3",
                "Ly6c2hi Mono" = "#27b5bf",
                "cDC2" = "#D3D3D3",
                "T-Cells.NK Cells" = "#D3D3D3",
                "IFN Mac" = "#d76da8")









## Subset
## ---------------------------------- #
R <- subset(R, cell_type %in% include_celltype )
root_cell <- "Ly6c2hi Mono"
Idents(R) <- R$cell_type %>% as.character()










## colors
## ---------------------------------- #
palette = plasma(100)
palette_celltype = brewer.pal(n = 8, name = "Dark2")

C <- palette_celltype[as.factor(R$cell_type)]
names(C) <- as.factor(R$cell_type)

color <- unique(C)
names(color) <- unique(names(C))








## run Slingshots 
## ---------------------------------- #
start.clus <- root_cell
reduction = 'umap'
sds = slingshot(Embeddings(R, reduction), clusterLabels = Idents(R), start.clus = start.clus )
R@tools[['slingshot']] = SlingshotDataSet(sds)
pseudotime = slingPseudotime(sds)


sds_all = slingshot(Embeddings(R_ALL, reduction), clusterLabels = Idents(R_ALL), start.clus = start.clus )








## Plot slingshot curves
## ---------------------------------- # 
pseudotime %>% as.data.frame() %>%
                        write.table(paste0(outdir, "slingshot.pseudotime.txt"), sep="\t", row.names = TRUE, col.names = TRUE)
curves = colnames(pseudotime)













## Plot slingshot curves
## ---------------------------------- #
sds_all$reducedDim %>% rownames() -> cell_id
umap_colors[R_ALL$cell_type] -> cell_color
names(cell_color) <- colnames(R_ALL)
cell_color[cell_id] -> U_C


pdf(file = paste0(outdir, 'slingshot_curves.all.umap.NEW.pdf'), width = 7.2)
print(plot(sds_all$reducedDim, col = U_C, pch = 16, ccex = 0.5) +
                lines(SlingshotDataSet(sds), lwd = 2, col = 'black')
                )
dev.off()










## Plot slingshot curve II by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds




R_TIME <- pseudotime_orig %>% as.data.frame() %>% dplyr::mutate("cell_id" = rownames(.))
R_META <- R_ALL[["cell_type"]] %>% dplyr::mutate("cell_id" = rownames(.))
R_UMAP <- Embeddings(object = R_ALL, reduction = "umap") %>% data.frame() %>% 
				dplyr::mutate("cell_id" = rownames(.)) %>% 
				left_join(., R_META, by="cell_id") %>% 
				left_join(., R_TIME, by="cell_id") 








make_umap <- function(df, color_column, color_list) {
	pdf(file = paste0(outdir, 'slingshot_curves.separate.colored.by.', color_column, '.NEW.pdf'), width = 7)
	print(ggplot(df, aes(x=UMAP_1, y=UMAP_2, group = eval(parse(text=color_column)), colour = eval(parse(text=color_column)))) +
	          geom_point(size=0.8, alpha=0.7) + 
	          theme_light() + theme_classic() +
		  scale_colour_manual(values=color_list) + 
	          theme(legend.position = "none", panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
	          scale_alpha(guide = 'none'))# to remove extra legend 
	dev.off()
	}

make_umap(R_UMAP, "cell_type", umap_colors)









make_umap_c <- function(df, color_column) {
	pdf(file = paste0(outdir, 'slingshot_curves.separate.colored.by.', color_column, '.NEW.pdf'), width = 7)
        print(ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, colour = eval(parse(text=color_column)) )) +
                  geom_point(size=0.8, alpha=0.7) +
                  theme_light() + theme_classic() +
		  scale_colour_viridis_c(na.value="#D3D3D3", option = "C") +
                  theme(legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                  scale_alpha(guide = 'none'))# to remove extra legend 
        dev.off()
	}  

make_umap_c(R_UMAP, "Lineage1")
make_umap_c(R_UMAP, "Lineage2")














## Plot slingshot curve by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds
pdf(file = paste0(outdir, 'slingshot_curves.separate.colored.by.time.all.umap.NEW.pdf'), width = 14, height=7)
par(mfrow = c(1, 2))
for ( c_num in seq(1, length(curves))) {

                        sds <- sds[,curves[c_num]]
			
			sds$reducedDim %>% rownames() -> cell_id
			sds_all$reducedDim %>% rownames() %>% as.data.frame() -> all_cell_id
			colnames(all_cell_id) <- "cell_id"
				
				
			pseudotime = slingPseudotime(sds)
			colors = palette[cut(pseudotime[,1], breaks = 100)]
			names(colors) <- cell_id
			colors %>% as.data.frame() %>% rownames_to_column() -> colors_df	
			colnames(colors_df) <- c("cell_id", "color")
			
			all_cell_id %>% 
				left_join(., colors_df, by="cell_id") %>% 
				mutate(color = ifelse(is.na(.[["color"]])==TRUE, "#D3D3D3", color)) -> new_colors 

                        print(plot(sds_all$reducedDim, col = new_colors$color, pch = 16, cex = 0.5, main = curves[c_num] ) +
                                lines(SlingshotDataSet(sds), linInd = c_num, lwd = 2, col = 'black'))


                        sds <- sds_orig
                        pseudotime <- pseudotime_orig



			}

dev.off()














## Add pseudotimes to meta data of R object 
## ---------------------------------- #
## > R[[c("Lineage1","Lineage2")]] %>% as.data.frame() -> b
## > pseudotime %>% as.data.frame() -> a
## > all.equal(a,b)
## [1] TRUE
## ---------------------------------- #
for ( curve in curves ) {
	pseudotime_sub <- pseudotime[colnames(R),curve]
	R <- AddMetaData(object = R,
                         metadata = pseudotime_sub,
                         col.name = curve
                         )
         }












## Condition density along pseudotime
## ---------------------------------- #
make_density <- function(your_obj, curve) {
		pdf(file = paste0(outdir, 'density_condition.', curve, '.pdf'), width = 7, height = 5)
		df <- data.frame(your_obj[["cell_type"]], your_obj[[curve]]) 
		colnames(df) <- c("cell_type", "Lineage")
		na.omit(df) -> df
		p <- ggplot(df, aes(x=Lineage, fill=cell_type)) +
			geom_density(alpha=0.4) + theme_classic()+
			scale_fill_manual(values=C) 
		print(p)
		dev.off()
		}

make_density(R, "Lineage1")
make_density(R, "Lineage2")












run_monocle_DE <- function(your_obj, your_column, VGAM_opt, VGAM_opt_name) {


	# subset based on lineage 2 
	# ---------------------------------- #
	L <- your_obj[[your_column]] %>% deframe()
	names(L) <- rownames(your_obj[[your_column]])
	L[!is.na(L)] %>% names() -> L2_cell
	#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
	your_obj[, L2_cell] -> new_R 
	Idents(new_R) <- new_R$cell_type %>% as.character()

	
	
		
	# Transfer into cds
	# ---------------------------------- #
	cds <- as.CellDataSet(new_R)
	# Estimate size factor
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	
	
	
	# call Monocle2
	# ---------------------------------- #
	# install https://www.bioconductor.org/packages/3.14/bioc/src/contrib/Archive/monocle/
	# refer to : http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
	# refer to : https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/monocle2.html#monocle2-process
	# ---------------------------------- #
	# select superset of feature genes as genes expressed in at least 5% of all the cells.
	# ---------------------------------- #
	cds <- detectGenes(cds, min_expr = 0.1)
	fData(cds)$use_for_ordering <- fData(cds)$num_cells_expressed > 0.05 * ncol(cds)
	cds_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))
	
	
	
	# get genes used for ordering cells 
	# ---------------------------------- #
	# while removing batch by using fullModelFormulaStr 
	# https://www.biostars.org/p/316204/
	# ---------------------------------- #
	clustering_DEG_genes <- differentialGeneTest(cds[cds_genes,],
						fullModelFormulaStr = '~orig.ident',
					     	cores = 10)





	cds_ordering_df <- clustering_DEG_genes %>% filter(qval < 0.01 & use_for_ordering == TRUE) %>% arrange(qval)
	cds_ordering_df[1:1000, ] %>% select(gene_short_name) %>% deframe() %>% as.character() -> cds_ordering_genes


	# Clustering Genes by Pseudotemporal Expression Pattern by Monocle2
	# ---------------------------------- #
        # df = 1 : a linear fit 
        # df = 2 : affords a little nonlinearity
        # df = 3 : VGAM
        # http://www2.uaem.mx/r-mirror/web/packages/VGAM/vignettes/categoricalVGAM.pdf
        # ---------------------------------- #  

	pData(cds)[[your_column]] -> pData(cds)$Pseudotime
	diff_test_res <- differentialGeneTest(cds[cds_ordering_genes,],
	                                fullModelFormulaStr = VGAM_opt, 
	                                cores = 10)
	
	diff_test_res %>% filter(qval <0.01) %>% arrange(qval) %>% write.table(., paste0(outdir, "slingshot.DEgenes.", your_column, ".", VGAM_opt_name, ".txt"), quote = FALSE, sep="\t", row.names = TRUE, col.names = TRUE)


	

	# make plots
	# ---------------------------------- #
	DE_N_Set <- c(50) # c(20, 30, 50)

	for ( DE_N in DE_N_Set ) { 
	
		diff_test_res %>% filter(qval <0.01) %>% arrange(qval) %>% top_n(., DE_N, wt=-qval) %>% rownames() -> sig_gene_names


	        hm <- get_pseudotime_matrix(cds[sig_gene_names,],  
	                                    cluster_rows = TRUE,
	                                    hclust_method = "ward.D",
	                                    num_clusters = 6,
	                                    hmcols = NULL,
	                                    add_annotation_row = NULL,
        	                            add_annotation_col = NULL,
        	                            show_rownames = FALSE,
        	                            use_gene_short_name = TRUE,
        	                            norm_method = "log",
        	                            scale_max=3,
        	                            scale_min=-3,
        	                            trend_formula = VGAM_opt, 
        	                            return_heatmap=TRUE,
        	                            cores=1)


		bks = c(seq(min(hm), 0, length.out=ceiling(200/2) + 1),
	                seq(max(hm)/200, max(hm),length.out=floor(200/2)))
		
		my_color4 = plasma(length(bks))
		my_color5 = colorRampPalette(rev(rcartocolor::carto_pal(7, "Sunset")))(length(bks))
		my_color6 = colorRampPalette(rcartocolor::carto_pal(7, "ag_Sunset"))(length(bks))
		my_color7 = colorRampPalette(rev(rcartocolor::carto_pal(7, "SunsetDark")))(length(bks))
			
		my_color_set <- list(my_color4, my_color5, my_color6, my_color7)
		my_color_name <- c("plasma", "Sunset", "ag_Sunset", "SunsetDark")
		


		# cluster and re-order rows
		# ---------------------------------- #
		# ALL_HCS <- c( "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
		# ---------------------------------- #

		ALL_HCS <- c("ward.D")

		for ( sub_color in seq(1,length(my_color_set)))  {
        		for ( ALL_HC in c(ALL_HCS) ) {	
		
				pdf(file = paste0(outdir, 'slingshot.top', DE_N , '.DEgenes.', your_column, '.', ALL_HC, '.', my_color_name[sub_color], '.', VGAM_opt_name, '.pdf'), width = 4, height = 6)
				print(monocle::plot_pseudotime_heatmap(cds[sig_gene_names,],
	        	        	        ##add_annotation_col = your_column,
						cluster_rows = TRUE,
						trend_formula = VGAM_opt,
						hclust_method = ALL_HC, 
						num_clusters = 1,
						hmcols = my_color_set[sub_color][[1]],
						scale_max = 3, 
						scale_min = -3,
	        	        	        cores = 1,
	        	        	        show_rownames = T,
						return_heatmap = FALSE))
				dev.off()
				}


			}

		}
	

	
	
	

        colors = palette[cut(pData(cds)$Pseudotime, breaks = 100)]
	phenoData(cds)[["color"]] <- colors
	GENE_OF_INTEREST <- c("Spp1", "Arg1", "Pf4", "Fn1")
	pdf(file = paste0(outdir, 'slingshot.gene_of_interest.pseudotime.', your_column,  '.plasma.pdf'), width = 4, height = 10)
	print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = your_column ) +         
	      				scale_color_viridis(option = "C") 
	      				#scale_color_viridis()
	      				)
	dev.off()

	pdf(file = paste0(outdir, 'slingshot.gene_of_interest.pseudotime.', your_column,  '.inferno.pdf'), width = 4, height = 10)
        print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = your_column ) +
                                        scale_color_viridis(option = "B")
                                        #scale_color_viridis()
                                        )

	dev.off()


	}



##run_monocle_DE(R, "Lineage1", "~sm.ns(Pseudotime, df=3)", "vgam")
##run_monocle_DE(R, "Lineage1", "~sm.ns(Pseudotime, df=2)", "non_linear")
run_monocle_DE(R, "Lineage1", "~sm.ns(Pseudotime, df=1)", "linear")

