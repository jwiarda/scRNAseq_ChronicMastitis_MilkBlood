library(Seurat)
library(future)
library(SeuratDisk)
library(ggplot)
library(writexl)
library(readxl)

set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly.h5seurat')
seu
DefaultAssay(seu) <- 'RNA'

# Cluster vs all (granulocytes only) ----
Idents(seu) <- seu$phyloorder

de <- FindAllMarkers(seu)
de <- subset(de, p_val_adj < 0.05) # make sure the adjusted p-values are still < 0.05 since some genes in DE list have p_val_adj > 0.05

# Save DE gene list:
write_xlsx(x = de, 
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/FinalClustersDEGs_GranulocytesOnly.xlsx",
           col_names = TRUE)

# Save background gene list:
bg <- data.frame(rownames(seu))
colnames(bg) <- 'gene'
write_xlsx(x = bg, 
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/BackgroundGeneList_GranulocytesOnly.xlsx",
           col_names = TRUE)

# Make heatmap of top DEGs:
#de <- read_excel('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/FinalClustersDEGs.xlsx')
de <- subset(de, avg_log2FC > 0) # only take genes enriched in the clusters 
Idents(seu) <- seu$phyloorder
topgenes <- de %>% group_by(cluster) %>% top_n(5, avg_log2FC) # only plot top 5 genes per cluster, as determined by highest average logFC values
DoHeatmap(subset(seu, downsample = 100), # take only 50 cells per cluster for plotting
          features = as.character(topgenes$gene), 
          assay = "RNA", 
          disp.min = -1, 
          disp.max = 2) +
  scale_fill_gradientn(colors = c("darkturquoise", "grey90", "indianred1", "red"))

# Perform DGE analysis between blood and milk within clusters ----


# Perform pairwise DGE for only granulocyte clusters ----

Idents(seu) <- seu$celltype
seu <- subset(seu, idents = 'granulocyte')
Idents(seu) <- seu$phyloorder
clusters <- unique(Idents(seu)) # identify all of our cluster IDs
pairwise <- combn(clusters, 2) # create all pairwise combinations of cluster IDs
p1 <- pairwise[1,] 
p2 <- pairwise[2,] 
comps1 <- data.frame(p1, p2)
colnames(comps1) <- c('pop1', 'pop2')
comps2 <- data.frame(p2, p1)
colnames(comps2) <- c('pop1', 'pop2')
comps <- rbind(comps1, comps2)

results <- list()
for(i in 1:nrow(comps)) {
  markers <- FindMarkers(seu, 
                         ident.1 = comps[i,1], 
                         ident.2 = comps[i,2],
                         assay = "RNA",
                         only.pos = TRUE) # only take positive logFC results since we perform recipricol comparisons for each cluster combination
  markers$gene <- rownames(markers)
  markers$pop1 <- paste(comps[i,1])
  markers$pop2 <- paste(comps[i,2])
  markers$comparison <- paste(markers$pop1, markers$pop2, sep = 'v')
  results[[i]] <- markers
} # if any of the comparisons don't turn up DE genes, this function won't work... it's then also likely the data has been over-clustered in preceding steps....
pwAll <- do.call(rbind, results)
pwAll <- subset(pwAll, p_val_adj < 0.05)

# Save DE gene list:
write_xlsx(x = pwAll, 
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/GranulocyteDEGs_pairwise.xlsx",
           col_names = TRUE)

up <- data.frame(table(pwAll$gene, pwAll$pop1)) # how many times does each gene show up as sig up in each cluster?
up <- subset(up, Freq == (length(clusters)-1)) # these are the genes that were highly specific to particular clusters since they are sig increased relative to every other pairwise cluster comparison

down <- data.frame(table(pwAll$gene, pwAll$pop2)) # how many times does each gene show up as sig down in each cluster?
down <- subset(down, Freq == (length(clusters)-1)) # these are the genes that were highly specific to particular clusters since they are sig decreased relative to every other pairwise cluster comparison

up30 <- subset(up, Var2 == 'c30')
up60 <- subset(up, Var2 == 'c60')
up46 <- subset(up, Var2 == 'c46')
up19 <- subset(up, Var2 == 'c19')

down30 <- subset(down, Var2 == 'c30')
down60 <- subset(down, Var2 == 'c60')
down46 <- subset(down, Var2 == 'c46')
down19 <- subset(down, Var2 == 'c19')

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] readxl_1.4.3          writexl_1.4.2         SeuratDisk_0.0.0.9020 future_1.33.1         SeuratObject_4.1.3    Seurat_4.3.0.1       

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-9           ellipsis_0.3.2         ggridges_0.5.4         rstudioapi_0.15.0      spatstat.data_3.0-1   
#[8] leiden_0.4.3           listenv_0.9.1          farver_2.1.1           ggrepel_0.9.5          bit64_4.0.5            fansi_1.0.6            codetools_0.2-19      
#[15] splines_4.2.2          polyclip_1.10-4        jsonlite_1.8.8         ica_1.0-3              cluster_2.1.4          png_0.1-8              uwot_0.1.16           
##[22] shiny_1.8.0            sctransform_0.3.5      spatstat.sparse_3.0-2  compiler_4.2.2         httr_1.4.7             Matrix_1.6-1           fastmap_1.1.1         
#[29] lazyeval_0.2.2         limma_3.52.4           cli_3.6.2              later_1.3.2            htmltools_0.5.7        tools_4.2.2            igraph_1.5.1          
#[36] gtable_0.3.4           glue_1.7.0             RANN_2.6.1             reshape2_1.4.4         dplyr_1.1.4            Rcpp_1.0.12            scattermore_1.2       
#[43] cellranger_1.1.0       vctrs_0.6.5            spatstat.explore_3.2-3 nlme_3.1-163           progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.1-6 
#[50] stringr_1.5.1          globals_0.16.2         mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3         
#[57] MASS_7.3-60            zoo_1.8-12             scales_1.3.0           promises_1.2.1         spatstat.utils_3.0-3   parallel_4.2.2         RColorBrewer_1.1-3    
#[64] reticulate_1.35.0      pbapply_1.7-2          gridExtra_2.3          ggplot2_3.5.0          stringi_1.8.3          rlang_1.1.3            pkgconfig_2.0.3       
#[71] matrixStats_1.0.0      lattice_0.21-8         ROCR_1.0-11            purrr_1.0.2            tensor_1.5             patchwork_1.2.0        htmlwidgets_1.6.4     
#[78] labeling_0.4.3         cowplot_1.1.3          bit_4.0.5              tidyselect_1.2.0       parallelly_1.37.1      RcppAnnoy_0.0.21       plyr_1.8.9            
#[85] magrittr_2.0.3         R6_2.5.1               generics_0.1.3         pillar_1.9.0           withr_3.0.0            fitdistrplus_1.1-11    survival_3.5-7        
#[92] abind_1.4-5            sp_2.0-0               tibble_3.2.1           future.apply_1.11.1    crayon_1.5.2           hdf5r_1.3.8            KernSmooth_2.23-22    
#[99] utf8_1.2.4             spatstat.geom_3.2-5    plotly_4.10.4          grid_4.2.2             data.table_1.15.2      digest_0.6.34          xtable_1.8-4          
#[106] tidyr_1.3.1            httpuv_1.6.14          munsell_0.5.0          viridisLite_0.4.2     