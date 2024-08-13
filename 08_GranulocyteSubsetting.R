library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230130_JEW_IntegratedSeurat_AnnotatedClusters_PhyloOrder.h5seurat')
DefaultAssay(seu) <- 'RNA'
Idents(seu) <- seu$celltype
seu <- subset(seu, idents = 'granulocyte')

DefaultAssay(seu) <- 'RNA'
counts <- as.data.frame(seu[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])
seu <- DietSeurat(seu, 
                  counts = TRUE,
                  data = TRUE,
                  scale.data = FALSE, # remove the scaled data
                  dimreducs = c('tsne', 'umap'),
                  features = keep, # keep only genes with non-zero counts across all cells
                  assays = 'RNA') # keep only RNA assay and remove SCT and integrated

seu.list <- SplitObject(seu, split.by = "sample_ID") # split by sample IDs

## Perform SCTransform normalization of data from each sample:
for (i in 1:length(seu.list)) { # normalize data using SCTransform method
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

## Integrate the data from different samples:
seu.features <- SelectIntegrationFeatures(seu.list, # select the genes to use for integration
                                          verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list, 
                               anchor.features = seu.features,
                               verbose = TRUE)
seu.anchors <- FindIntegrationAnchors(seu.list, # identify anchors for integration from top PCs
                                      normalization.method = "SCT", 
                                      anchor.features = seu.features)
seu.integrated <- IntegrateData(seu.anchors, # integrate data
                                normalization.method = "SCT")

## Run multidimensional analyses on data:
seu.integrated <- RunPCA(seu.integrated, # run PCA analysis for 100 dimensions of the data
                         npcs = 100, 
                         verbose = TRUE) 
ElbowPlot(seu.integrated,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering


## Also define nearest neighbors:
seu.integrated <- FindNeighbors(seu.integrated, 
                                dims = 1:30)

 ## Add normalized/scaled data to RNA assay:
#dim(seu.integrated[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
seu.integrated <- NormalizeData(seu.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
seu.integrated <- ScaleData(seu.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA")
#dim(seu.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now

# Reassign UMAP/t-SNE coordinates
identical(colnames(seu), colnames(seu.integrated)) # make sure all cell barcodes line up first
seu.integrated@reductions$umap <- seu@reductions$umap
seu.integrated@reductions$tsne <- seu@reductions$tsne
                                
DefaultAssay(seu.integrated) <- "RNA"
SaveH5Seurat(seu.integrated, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly.h5seurat', overwrite = TRUE)
     
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
#[22] shiny_1.8.0            sctransform_0.3.5      spatstat.sparse_3.0-2  compiler_4.2.2         httr_1.4.7             Matrix_1.6-1           fastmap_1.1.1         
#[29] lazyeval_0.2.2         cli_3.6.2              later_1.3.2            htmltools_0.5.7        tools_4.2.2            igraph_1.5.1           gtable_0.3.4          
#[36] glue_1.7.0             RANN_2.6.1             reshape2_1.4.4         dplyr_1.1.4            Rcpp_1.0.12            scattermore_1.2        cellranger_1.1.0      
#[43] vctrs_0.6.5            spatstat.explore_3.2-3 nlme_3.1-163           progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.1-6  stringr_1.5.1         
#[50] globals_0.16.2         mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3          MASS_7.3-60           
#[57] zoo_1.8-12             scales_1.3.0           promises_1.2.1         spatstat.utils_3.0-3   parallel_4.2.2         RColorBrewer_1.1-3     reticulate_1.35.0     
#[64] pbapply_1.7-2          gridExtra_2.3          ggplot2_3.5.0          stringi_1.8.3          rlang_1.1.3            pkgconfig_2.0.3        matrixStats_1.0.0     
#[71] lattice_0.21-8         ROCR_1.0-11            purrr_1.0.2            tensor_1.5             patchwork_1.2.0        htmlwidgets_1.6.4      labeling_0.4.3        
#[78] cowplot_1.1.3          bit_4.0.5              tidyselect_1.2.0       parallelly_1.37.1      RcppAnnoy_0.0.21       plyr_1.8.9             magrittr_2.0.3        
#[85] R6_2.5.1               generics_0.1.3         pillar_1.9.0           withr_3.0.0            fitdistrplus_1.1-11    survival_3.5-7         abind_1.4-5           
#[92] sp_2.0-0               tibble_3.2.1           future.apply_1.11.1    crayon_1.5.2           hdf5r_1.3.8            KernSmooth_2.23-22     utf8_1.2.4            
#[99] spatstat.geom_3.2-5    plotly_4.10.4          grid_4.2.2             data.table_1.15.2      digest_0.6.34          xtable_1.8-4           tidyr_1.3.1           
#[106] httpuv_1.6.14          munsell_0.5.0          viridisLite_0.4.2   