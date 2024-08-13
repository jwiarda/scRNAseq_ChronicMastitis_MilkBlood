library(Seurat)
library(future)
library(SeuratDisk)
library(clustree)
library(ggplot)
library(dplyr)

set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230127_JEW_IntegratedSeurat.h5seurat')
seu

DefaultAssay(seu) <- 'integrated'
seu <- FindClusters(seu, 
                    resolution = c(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10), 
                    verbose = FALSE) 

clustree(seu, # observe patterns of clustering with clustree
         prefix = "integrated_snn_res.")
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/clustree.jpeg')

DefaultAssay(seu) <- 'RNA'
FeaturePlot(seu,
            features = c('PTPRC', # leukocyte
                         'LALBA', 'CSN2', # non-immune
                         'CD3E', 'CD247', 'ZAP70', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', # T/ILC
                         'CD79A', 'CD19', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 'XBP1', # B/ASC
                         'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                         'CSF1R', 'CST3', 'CD83', 'CD86', 'CD14', 'TLR4', # macrophage/monocyte
                         'FCGR3A', # CD16 monocytes
                         'FLT3', 'ENSBTAG00000009656', 'ENSBTAG00000013919', # DC
                         'GZMB', 'IL3RA', # pDC...also CD4+ CD14- 
                         'CSF3R', 'CXCL8', 'SRGN', #granulocyte
                         'PCLAF', 'UBE2C', 'MKI67', # cell cycle
                         'AHSP', 'HBM'), # RBC genes
            cols = c('grey90', 'darkslateblue'),
            ncol = 11, 
            #order = TRUE, # pulls cells expressing genes to front
            reduction = 'tsne') & NoLegend() & NoAxes()
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/FeaturePlot_AllCells_CanonicalGenes.jpeg')

#DimPlot(seu, group.by = 'integrated_snn_res.4', reduction = 'tsne', label = TRUE)

DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
            features = c('PTPRC', # leukocyte
                         'LALBA', 'CSN2', # non-immune
                         'CD3E', 'CD247', 'ZAP70', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', # T/ILC
                         'CD79A', 'CD19', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 'XBP1', # B/ASC
                         'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                         'CSF1R', 'CST3', 'CD83', 'CD86', 'CD14', 'TLR4', # macrophage/monocyte
                         'FCGR3A', # CD16 monocytes
                         'FLT3', 'ENSBTAG00000009656', 'ENSBTAG00000013919', # DC
                         'GZMB', 'IL3RA', # pDC...also CD4+ CD14- 
                         'CSF3R', 'CXCL8', 'SRGN', #granulocyte
                         'PCLAF', 'UBE2C', 'MKI67', # cell cycle
                         'AHSP', 'HBM'), # RBC genes
            cols = c('goldenrod1', 'firebrick'),
        group.by = 'integrated_snn_res.4') & RotatedAxis()
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DotPlot_AllCells_CanonicalGenes_res4Clusters.jpeg')

seu$celltype <- seu$integrated_snn_res.4
Idents(seu) <- seu$celltype
types <- c(rep('granulocyte', 12), 'B/ASC', rep('T/ILC', 2), 'granulocyte',
           'T/ILC', rep('granulocyte', 2), 'monocyte/macrophage/cDC', 'B/ASC',
           'T/ILC', rep('granulocyte', 4), 'B/ASC', rep('granulocyte', 3),
           'B/ASC', rep('granulocyte', 2), rep('monocyte/macrophage/cDC', 2),
           'B/ASC', 'T/ILC', 'monocyte/macrophage/cDC', 'granulocyte', 'T/ILC', 
           'B/ASC', 'T/ILC', 'monocyte/macrophage/cDC', 'non-immune', 'T/ILC', 'granulocyte',
           'T/ILC', 'granulocyte', 'T/ILC', 'B/ASC', 'monocyte/macrophage/cDC',
           'granulocyte', 'pDC', 'cycling B/myeloid', 'monocyte/macrophage/cDC', 'T/ILC', 
           'monocyte/macrophage/cDC', 'granulocyte', 'monocyte/macrophage/cDC', 'granulocyte', 
           'T/ILC') # Rename clusters based on phenotype IDs
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$celltype <- Idents(seu)
Idents(seu) <- seu$celltype
#DimPlot(seu, reduction = 'tsne', label = TRUE)

seu2 <- subset(seu, idents = c('cycling B/myeloid'))
DefaultAssay(seu2) <- 'integrated'
seu2 <- RunPCA(seu2, # run PCA analysis for 100 dimensions of the data
                         npcs = 100, 
                         verbose = TRUE) 
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2, 
                    resolution = c(.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10), 
                    verbose = FALSE) 

DefaultAssay(seu2) <- 'RNA'
FeaturePlot(seu2,
            features = c(
                         'CD79A', 'PAX5', 'JCHAIN', # B/ASC
                         'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                         'CSF1R', 'CD14', 'TLR4'), # macrophage/monocyte
            cols = c('grey90', 'darkslateblue'),
            ncol = 3, 
            #order = TRUE, # pulls cells expressing genes to front
            reduction = 'tsne') & NoLegend() & NoAxes()
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/FeaturePlot_CycBMyeloid_CanonicalGenes.jpeg')

DefaultAssay(seu2) <- 'RNA'
#DimPlot(seu2, reduction = 'tsne', label = TRUE, group.by = 'integrated_snn_res.2')
DotPlot(seu2,
        features = c('CD79A', 'PAX5', 'JCHAIN', # B/ASC
                     'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                     'CSF1R', 'CD14', 'TLR4'), # macrophage/monocyte
        cols = c('goldenrod1', 'firebrick'),
        group.by = 'integrated_snn_res.2') & RotatedAxis()
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DotPlot_CycBMyeloid_CanonicalGenes_res2Clusters.jpeg')


seu2$celltype <- seu2$integrated_snn_res.2
Idents(seu2) <- seu2$celltype
types <- c(rep('B/ASC', 2), 'monocyte/macrophage/cDC', rep('B/ASC', 3), 
           rep('monocyte/macrophage/cDC', 2)) # Rename clusters based on phenotype IDs
names(types) <- levels(seu2) # assign GutCellTypes to cluster numbers
seu2 <- RenameIdents(seu2, types) # change dataset identity to cell types in new Seurat object
seu2$celltype <- Idents(seu2)
Idents(seu2) <- seu2$celltype
#DimPlot(seu2, reduction = 'tsne')

DotPlot(seu2,
        features = c('CD79A', 'PAX5', 'JCHAIN', # B/ASC
                     'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                     'CSF1R', 'CD14', 'TLR4'), # macrophage/monocyte
        cols = c('goldenrod1', 'firebrick'),
        group.by = 'celltype') & RotatedAxis()
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DotPlot_CycBMyeloid_CanonicalGenes_CalledCells.jpeg')


SaveH5Seurat(seu2, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230127_JEW_CycBMyeloid_sub.h5seurat', overwrite = TRUE)

## Incorporate subset cluster 58 annotations to make new cluster IDs at clustering resolution of 4:

cycM <- rownames(seu2@meta.data %>% filter(seu2$celltype == 'monocyte/macrophage/cDC'))
cycB <- rownames(seu2@meta.data %>% filter(seu2$celltype == 'B/ASC'))

Idents(seu) <- seu$integrated_snn_res.4
Idents(seu, cells = cycM) <- 'cycM'
Idents(seu, cells = cycB) <- 'cycB'
seu$cluster <- Idents(seu)

seu$celltype <- seu$cluster
Idents(seu) <- seu$celltype
types <- c('B/ASC', 'monocyte/macrophage/cDC', rep('granulocyte', 12), 'B/ASC', rep('T/ILC', 2), 'granulocyte',
           'T/ILC', rep('granulocyte', 2), 'monocyte/macrophage/cDC', 'B/ASC',
           'T/ILC', rep('granulocyte', 4), 'B/ASC', rep('granulocyte', 3),
           'B/ASC', rep('granulocyte', 2), rep('monocyte/macrophage/cDC', 2),
           'B/ASC', 'T/ILC', 'monocyte/macrophage/cDC', 'granulocyte', 'T/ILC', 
           'B/ASC', 'T/ILC', 'monocyte/macrophage/cDC', 'non-immune', 'T/ILC', 'granulocyte',
           'T/ILC', 'granulocyte', 'T/ILC', 'B/ASC', 'monocyte/macrophage/cDC',
           'granulocyte', 'pDC', 'monocyte/macrophage/cDC', 'T/ILC', 
           'monocyte/macrophage/cDC', 'granulocyte', 'monocyte/macrophage/cDC', 'granulocyte', 
           'T/ILC') # Rename clusters based on phenotype IDs
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$celltype <- Idents(seu)
Idents(seu) <- seu$celltype

seu$cluster <- seu$cluster
Idents(seu) <- seu$cluster
types <- c('c54', 'c62', paste0('c', 1:53), paste0('c', 55:61)) # Rename clusters based on phenotype IDs
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$cluster <- Idents(seu)
Idents(seu) <- seu$cluster
seu$cluster <- factor(seu$cluster,levels=c(paste0('c', 1:62)))

Idents(seu) <- seu$cluster
DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = c('PTPRC', # leukocyte
                     'LALBA', 'CSN2', # non-immune
                     'CD3E', 'CD247', 'ZAP70', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', # T/ILC
                     'CD79A', 'CD19', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 'XBP1', # B/ASC
                     'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                     'CSF1R', 'CST3', 'CD83', 'CD86', 'CD14', 'TLR4', # macrophage/monocyte
                     'FCGR3A', # CD16 monocytes
                     'FLT3', 'ENSBTAG00000009656', 'ENSBTAG00000013919', # DC
                     'GZMB', 'IL3RA', # pDC...also CD4+ CD14- MHCII-
                     'CSF3R', 'CXCL8', 'SRGN', #granulocyte
                     'PCLAF', 'UBE2C', 'MKI67', # cell cycle
                     'AHSP', 'HBM'), # RBC genes
        cols = c('goldenrod1', 'firebrick'),
        group.by = 'cluster') & RotatedAxis()
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DotPlot_AllCells_CanonicalGenes_FinalClusters.jpeg')

DimPlot(seu, group.by = 'sample_ID', reduction = 'tsne', shuffle = TRUE) #& NoLegend()
#ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_SampleID_NoLegend.jpeg')
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_SampleID_Legend.jpeg')

DimPlot(seu, group.by = 'tissue', reduction = 'tsne', shuffle = TRUE, cols = c('indianred1', 'royalblue')) #& NoLegend()
#ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_tissue_NoLegend.jpeg')
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_tissue_Legend.jpeg')

DimPlot(seu, reduction = 'tsne', group.by = 'celltype', shuffle = TRUE,
        cols = c('slateblue', 'violetred', 'chartreuse4', 'chocolate2', 'goldenrod1', 'deepskyblue2', 'grey80')) #& NoLegend()
#ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_celltype_NoLegend.jpeg')
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_celltype_Legend.jpeg')

DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'integrated_snn_res.4', shuffle = TRUE, repel = TRUE) #& NoLegend()
#ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_res4clusters_NoLegend.jpeg')
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_res4clusters_Legend.jpeg')

DimPlot(seu, reduction = 'tsne', label = TRUE, group.by = 'cluster', shuffle = TRUE, repel = TRUE) #& NoLegend()
#ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_finalClusters_NoLegend.jpeg')
ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DimPlot_AllCells_finalClusters_Legend.jpeg')

FeaturePlot(seu, 
            cols = c('grey90', 'darkslateblue'),
            reduction = 'tsne',
            features = c('CD79A', 'CSF1R', 'CSF3R', 'CD3E', 'LALBA', 'CSN2', 'IL3RA')) & NoAxes()

FeaturePlot(seu,
            features = c('nCount_RNA', 'nFeature_RNA'), 
            max.cutoff = 10000, 
            reduction = 'tsne',
            cols = c('grey90', 'darkslateblue'))

SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230129_JEW_IntegratedSeurat_AnnotatedClusters.h5seurat', overwrite = TRUE)

# Jerry-rig data names to be compatible with clustree graphing:
seu$plot.1 <- seu$celltype
seu$plot.2 <- seu$cluster
seu$plot.3 <- seu$integrated_snn_res.4
clustree(seu, # observe patterns of clustering with clustree
         prefix = 'plot.')

# Final dot plot:

Idents(seu) <- seu$cluster
levels(seu) <- c('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 
                 'c12', 'c16', 'c18', 'c19', 'c23', 'c24', 'c25', 'c26', 
                 'c28', 'c29', 'c30', 'c32', 'c33', 'c39', 'c46', 'c48', 'c52', 'c58', 'c60',
                 
                 'c20', 'c34', 'c35', 'c38', 'c43', 'c51', 'c55', 'c57', 'c59', 'c62',
                 
                 'c13', 'c21', 'c27', 'c31', 'c36', 'c41', 'c50', 'c54',
                 
                 'c14', 'c15', 'c17', 'c22', 'c37', 'c40', 'c42', 'c45', 'c47', 'c49', 'c56', 'c61',
                 
                 'c53',
                 
                 'c44')

DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = c('SRGN', 'CSF3R', 'CXCL8', #granulocyte
                     'SIRPA', 'CD68', 'TYROBP', # pan-myeloid
                     'TLR4', 'CSF1R', 'CST3', 'CD86', 'CD14', 'CD83', # macrophage/monocyte
                     'FCGR3A', # CD16 monocytes
                     'FLT3', 'ENSBTAG00000009656', 'ENSBTAG00000013919', # DC
                     'CD79A', 'CD19', 'PAX5', 'JCHAIN', 'IRF4', 'PRDM1', 'XBP1', # B/ASC
                     'CD3E', 'CD247', 'ZAP70', 'NCR1', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', # T/ILC
                     'GZMB', 'IL3RA', # pDC...also CD4+ CD14- MHCII-
                     'PTPRC', # leukocyte
                     'LALBA', 'CSN2', # non-immune
                     'AHSP', 'HBM'), # RBC genes
        cols = c('goldenrod1', 'firebrick')) & RotatedAxis()

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.1 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] dplyr_1.0.10          clustree_0.5.0        ggraph_2.0.6          ggplot2_3.3.6         SeuratDisk_0.0.0.9020 future_1.28.0        
#[7] sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         

#loaded via a namespace (and not attached):
#  [1] systemfonts_1.0.4           plyr_1.8.7                  igraph_1.3.4                lazyeval_0.2.2              splines_4.2.2              
#[6] BiocParallel_1.30.3         listenv_0.8.0               scattermore_0.8             GenomeInfoDb_1.32.4         digest_0.6.29              
#[11] htmltools_0.5.3             viridis_0.6.2               fansi_1.0.3                 magrittr_2.0.3              ScaledMatrix_1.4.1         
#[16] tensor_1.5                  cluster_2.1.4               ROCR_1.0-11                 limma_3.52.3                globals_0.16.1             
#[21] graphlayouts_0.8.1          matrixStats_0.62.0          spatstat.sparse_2.1-1       colorspace_2.0-3            ggrepel_0.9.1              
#[26] textshaping_0.3.6           xfun_0.33                   crayon_1.5.1                RCurl_1.98-1.8              jsonlite_1.8.0             
#[31] spatstat.data_2.2-0         progressr_0.11.0            survival_3.4-0              zoo_1.8-10                  ape_5.6-2                  
#[36] glue_1.6.2                  polyclip_1.10-0             gtable_0.3.1                zlibbioc_1.42.0             XVector_0.36.0             
#[41] leiden_0.4.3                DelayedArray_0.22.0         BiocSingular_1.12.0         future.apply_1.9.1          SingleCellExperiment_1.18.0
#[46] BiocGenerics_0.42.0         abind_1.4-5                 scales_1.2.1                DBI_1.1.3                   edgeR_3.38.4               
#[51] spatstat.random_2.2-0       miniUI_0.1.1.1              Rcpp_1.0.9                  viridisLite_0.4.1           xtable_1.8-4               
#[56] spatstat.core_2.4-4         reticulate_1.26             bit_4.0.4                   rsvd_1.0.5                  stats4_4.2.2               
#[61] htmlwidgets_1.5.4           httr_1.4.4                  RColorBrewer_1.1-3          ellipsis_0.3.2              ica_1.0-3                  
#[66] pkgconfig_2.0.3             farver_2.1.1                uwot_0.1.14                 deldir_1.0-6                locfit_1.5-9.6             
#[71] utf8_1.2.2                  labeling_0.4.2              reshape2_1.4.4              tidyselect_1.1.2            rlang_1.0.6                
#[76] later_1.3.0                 munsell_0.5.0               tools_4.2.2                 cli_3.4.0                   generics_0.1.3             
#[81] ggridges_0.5.3              evaluate_0.16               stringr_1.4.1               fastmap_1.1.0               ragg_1.2.2                 
#[86] goftest_1.2-3               yaml_2.3.5                  bit64_4.0.5                 knitr_1.40                  fitdistrplus_1.1-8         
#[91] tidygraph_1.2.2             purrr_0.3.4                 RANN_2.6.1                  pbapply_1.5-0               nlme_3.1-159               
#[96] mime_0.12                   hdf5r_1.3.5                 compiler_4.2.2              rstudioapi_0.14             beeswarm_0.4.0             
#[101] plotly_4.10.0               png_0.1-7                   spatstat.utils_2.3-1        tibble_3.1.8                tweenr_2.0.2               
#[106] stringi_1.7.8               rgeos_0.5-9                 lattice_0.20-45             Matrix_1.5-1                vctrs_0.4.1                
#[111] pillar_1.8.1                lifecycle_1.0.2             spatstat.geom_2.4-0         lmtest_0.9-40               RcppAnnoy_0.0.19           
#[116] BiocNeighbors_1.14.0        data.table_1.14.2           cowplot_1.1.1               bitops_1.0-7                irlba_2.3.5                
#[121] httpuv_1.6.6                patchwork_1.1.2             GenomicRanges_1.48.0        R6_2.5.1                    promises_1.2.0.1           
#[126] KernSmooth_2.23-20          gridExtra_2.3               vipor_0.4.5                 IRanges_2.30.1              parallelly_1.32.1          
#[131] codetools_0.2-18            MASS_7.3-58.1               gtools_3.9.3                assertthat_0.2.1            SummarizedExperiment_1.26.1
#[136] withr_2.5.0                 sctransform_0.3.4           S4Vectors_0.34.0            GenomeInfoDbData_1.2.8      mgcv_1.8-40                
#[141] parallel_4.2.2              miloR_1.4.0                 rpart_4.1.16                grid_4.2.2                  beachmat_2.12.0            
#[146] tidyr_1.2.1                 rmarkdown_2.16              MatrixGenerics_1.8.1        Rtsne_0.16                  ggforce_0.3.4              
#[151] Biobase_2.56.0              shiny_1.7.2                 ggbeeswarm_0.6.0     