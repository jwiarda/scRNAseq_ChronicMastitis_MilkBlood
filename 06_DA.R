library(miloR)
library(SingleCellExperiment)
library(data.table)
library(ggplot2)
library(dplyr)
library(writexl)
library(Seurat)
library(SeuratDisk)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230130_JEW_IntegratedSeurat_AnnotatedClusters_PhyloOrder.h5seurat')
Idents(seu) <- seu$cluster
seu$cluster <- factor(seu$cluster, levels = c('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 
                                              'c12', 'c16', 'c18', 'c19', 'c23', 'c24', 'c25', 'c26', 
                                              'c28', 'c29', 'c30', 'c32', 'c33', 'c39', 'c46', 'c48', 'c52', 'c58', 'c60',
                                              
                                              'c20', 'c34', 'c35', 'c38', 'c43', 'c51', 'c55', 'c57', 'c59', 'c62',
                                              
                                              'c13', 'c21', 'c27', 'c31', 'c36', 'c41', 'c50', 'c54',
                                              
                                              'c14', 'c15', 'c17', 'c22', 'c37', 'c40', 'c42', 'c45', 'c47', 'c49', 'c56', 'c61',
                                              
                                              'c53',
                                              
                                              'c44'))

set.seed(111) # set a seed for reproducibility
milo <- as.SingleCellExperiment(seu, assay = 'RNA')
milo <- Milo(milo)

#Also incorporate the shared nearest neighbors (SNN) graph calculated in Seurat into the Milo object. We will use this to calculate cell neighborhoods.

#Next, find and add neighbors:
#seu <- FindNeighbors(seu, assay = 'integrated', k.param = 20, dims = 1:30) # calculate nearest neighbors if seu@graphs$integrated_snn does not exist yet
miloR::graph(milo) <- miloR::graph(buildFromAdjacency(seu@graphs$integrated_snn, k=20))


### Create & visualize cell neighborhoods:

#Start by creating cell neighborhoods. Parameters of prop, k, and d may be modified slightly depending on the dataset used. Higher proportions of cells (prop) will take longer to run, but may require up to a value of 0.2 for smaller datasets. We choose to set k and d parameters according to those used to calculate our SNN graph and 'significant' PCs in our previous Seurat analysis.

#Now generate the neighborhoods:

set.seed(111) # set a seed for reproducibility of neighborhood generation
milo <- makeNhoods(milo,
                   prop = 0.2, # sample 20% of cells...probably safe to lower as far as 0.05 for datasets with >30k cells...may consider using proportions up to 0.2 if that helps optimize neighborhood size distribution peak; higher prop will increase run time. We will keep prop at 0.2 since we have some rare cell populations we want to make sure are represented.
                   k = 20, # set to k = 20 because for Seurat FindNeighbors() we used default k.param = 20
                   d = 30, # set to PCs used to find neighbors in Seurat
                   refined = TRUE, # always use refined unless you use graph-based data batch correction, then consider either-or
                   refinement_scheme="graph") # use graph-based approach so bottlenecking step calcNhoodDistance() isn't required

#Now that we have calculated our cell neighborhoods, let's look at their sizes. Ideally, peak size should fall between 50-100 cells per neighborhood but may be less for extremely small datasets:
plotNhoodSizeHist(milo) # ideally have peak of distribution between 50 and 100...otherwise consider increasing k or prop...peak may be <50 for small datasets

#Now let's move on to look at these cell neighborhoods overlaid onto UMAP coordinates:

milo <- buildNhoodGraph(milo)
plotNhoodGraph(milo, layout = 'TSNE')

### Count cells in each neighborhood

#Now let's do a head count of which cells came from each of our samples within each of our detected cell neighborhoods:

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="sample_ID")
head(nhoodCounts(milo))

### Create experimental design

#Create a model of experimental design variables:

milo_design <- data.frame(colData(milo))[,c("sample_ID", "tissue")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$sample_ID
milo_design

### Perform DA testing

#Perform DA testing on each neighborhood:

da_results <- testNhoods(milo,
                         design = ~ tissue,
                         design.df = milo_design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         reduced.dim = 'PCA')
head(da_results)

#Make a histogram of p-values found across cell neighborhoods:

ggplot(da_results, aes(PValue)) + geom_histogram(bins=100)

#Make a volcano plot of DA. Each dot is one cell neighborhood:

ggplot(da_results, aes(logFC, -log10(FDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

#Overlay logFC scores onto cell neighborhood central coordinates on t-SNE & UMAP plots:

plotNhoodGraphDA(milo, da_results, layout="TSNE",alpha=0.05)

#And we can also look at all cell neighborhoods on a bee swarm plot:

plotDAbeeswarm(da_results, alpha = 0.05)

#The problem right now is that we know which cell neighborhoods are DA and can guess at which cell types these correspond to, but we don't know for certain. Therefore, it may be helpful to annotate our cell neighborhoods and then re-assess DA....

### Annotate cell neighborhoods

#Annotate cell neighborhoods by finding which neighborhoods of most cells belonging to a specific cell ID. 

#Start by calculating the percentage of cells belonging to each annotated cell type within each neighborhood. We will record the annotation that has the largest percentage in each neighborhood:

da_results <- annotateNhoods(milo, da_results, coldata_col = "cluster")
da_results <- annotateNhoods(milo, da_results, coldata_col = "celltype")
head(da_results)

#Create a histogram to look at the largest percentages for a single cell type within each cell neighborhood:

ggplot(da_results, aes(cluster_fraction)) + geom_histogram(bins=100)
ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=100)

#Based on this graph, we need to set a cut-off value for what percentage of cells in a neighborhood must share a single ID to be assigned as that cell type annotation. In this case, we will set our cutoff at 0.7.

#Based on this criteria, any cell neighborhood with >70% of cells belonging to a single annotation will be assigned to that identity. Any cell neighborhood with <70% of cells belonging to a single annotation will be catergorized as 'Mixed' cell neighborhoods:

da_results$cluster <- ifelse(da_results$cluster_fraction <= 0.9, "Mixed", da_results$cluster)
da_results$celltype <- ifelse(da_results$celltype_fraction <= 0.9, "Mixed", da_results$celltype)

da_results$cluster <- factor(da_results$cluster, levels = levels(seu$cluster))
da_results$celltype <- factor(da_results$celltype, 
                              levels = c('non-immune', 'B/ASC', 'T/ILC', 'pDC',
                                         'monocyte/macrophage/cDC', 'granulocyte',
                                         'Mixed'))
### Plot DA across annotated cell neighborhoods:

#Make a bee swarm plot:

plotDAbeeswarm(da_results, group.by = "cluster", alpha = 0.05)
plotDAbeeswarm(da_results, group.by = "celltype", alpha = 0.05)

### Further summarization of DA results:

#Let's further summarize only those neighborhoods with an annotated cell type. Start by subsetting only non-mixed neighborhoods, creating a column defining DA significance (FDR < 0.05), and creating a column indicating fold-change towards enrichment in one group vs another:

da_sum <- subset(da_results, celltype_fraction > 0.9)
da_sum <- da_sum %>% mutate(significant_celltype = case_when(FDR < 0.05 ~ 'Sig', FDR >= 0.05 ~ 'NonSig'))
da_sum <- da_sum %>% mutate(FC = case_when(logFC > 0 ~ 'milk', logFC < 0 ~ 'blood', logFC == 0 ~ 'Neutral'))
da_sum$result_celltype <- paste(da_sum$FC, da_sum$significant_celltype, sep = '_')
da_sum$result_celltype <- replace(da_sum$result_celltype, da_sum$result_celltype == 'blood_NonSig', 'NonSig')
da_sum$result_celltype <- replace(da_sum$result_celltype, da_sum$result_celltype == 'milk_NonSig', 'NonSig')
table(da_sum$result_celltype, da_sum$celltype) # see summary of results per cell type

da_sum <- subset(da_results, cluster_fraction > 0.9)
da_sum <- da_sum %>% mutate(significant_cluster = case_when(FDR < 0.05 ~ 'Sig', FDR >= 0.05 ~ 'NonSig'))
da_sum <- da_sum %>% mutate(FC = case_when(logFC > 0 ~ 'milk', logFC < 0 ~ 'blood', logFC == 0 ~ 'Neutral'))
da_sum$result_cluster <- paste(da_sum$FC, da_sum$significant_cluster, sep = '_')
da_sum$result_cluster <- replace(da_sum$result_cluster, da_sum$result_cluster == 'blood_NonSig', 'NonSig')
da_sum$result_cluster <- replace(da_sum$result_cluster, da_sum$result_cluster == 'milk_NonSig', 'NonSig')
table(da_sum$result_cluster, da_sum$cluster) # see summary of results per cell type

### Save DA results:

write_xlsx(da_results, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DA/DAResults_AllCells.xlsx')

#Save Milo object:

saveRDS(milo, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DA/20230131_milo_AllCells.rds')

### View session information

sessionInfo()

#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] SingleCellExperiment_1.18.1 SummarizedExperiment_1.26.1 Biobase_2.56.0              GenomicRanges_1.48.0       
#[5] GenomeInfoDb_1.32.4         IRanges_2.30.1              S4Vectors_0.34.0            BiocGenerics_0.42.0        
#[9] MatrixGenerics_1.8.1        matrixStats_1.0.0           miloR_1.4.0                 edgeR_3.38.4               
#[13] limma_3.52.4                scales_1.3.0                ggplot2_3.5.0               stringr_1.5.1              
#[17] dplyr_1.1.4                 data.table_1.15.2           readxl_1.4.3                writexl_1.4.2              
#[21] SeuratDisk_0.0.0.9020       future_1.33.1               SeuratObject_4.1.3          Seurat_4.3.0.1             

#loaded via a namespace (and not attached):
#  [1] plyr_1.8.9             igraph_1.5.1           lazyeval_0.2.2         sp_2.0-0               splines_4.2.2         
#[6] BiocParallel_1.30.4    listenv_0.9.1          scattermore_1.2        digest_0.6.34          htmltools_0.5.7       
#[11] viridis_0.6.4          fansi_1.0.6            magrittr_2.0.3         ScaledMatrix_1.4.1     tensor_1.5            
#[16] cluster_2.1.4          ROCR_1.0-11            globals_0.16.2         graphlayouts_1.0.0     spatstat.sparse_3.0-2 
#[21] colorspace_2.1-0       ggrepel_0.9.5          crayon_1.5.2           RCurl_1.98-1.12        jsonlite_1.8.8        
#[26] progressr_0.14.0       spatstat.data_3.0-1    ape_5.7-1              survival_3.5-7         zoo_1.8-12            
#[31] glue_1.7.0             polyclip_1.10-4        gtable_0.3.4           zlibbioc_1.42.0        XVector_0.36.0        
#[36] leiden_0.4.3           DelayedArray_0.22.0    BiocSingular_1.12.0    future.apply_1.11.1    abind_1.4-5           
#[41] DBI_1.2.2              spatstat.random_3.1-6  miniUI_0.1.1.1         Rcpp_1.0.12            viridisLite_0.4.2     
#[46] xtable_1.8-4           reticulate_1.35.0      rsvd_1.0.5             bit_4.0.5              htmlwidgets_1.6.4     
#[51] httr_1.4.7             RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3              pkgconfig_2.0.3       
#[56] farver_2.1.1           uwot_0.1.16            deldir_1.0-9           locfit_1.5-9.8         utf8_1.2.4            
#[61] tidyselect_1.2.0       rlang_1.1.3            reshape2_1.4.4         later_1.3.2            munsell_0.5.0         
#[66] cellranger_1.1.0       tools_4.2.2            cli_3.6.2              generics_0.1.3         ggridges_0.5.4        
#[71] fastmap_1.1.1          goftest_1.2-3          bit64_4.0.5            fitdistrplus_1.1-11    tidygraph_1.2.3       
#[76] purrr_1.0.2            RANN_2.6.1             ggraph_2.1.0           pbapply_1.7-2          nlme_3.1-163          
#[81] mime_0.12              hdf5r_1.3.8            compiler_4.2.2         rstudioapi_0.15.0      beeswarm_0.4.0        
#[86] plotly_4.10.4          png_0.1-8              spatstat.utils_3.0-3   tibble_3.2.1           tweenr_2.0.2          
#[91] stringi_1.8.3          lattice_0.21-8         Matrix_1.6-1           vctrs_0.6.5            pillar_1.9.0          
#[96] lifecycle_1.0.4        spatstat.geom_3.2-5    lmtest_0.9-40          RcppAnnoy_0.0.21       BiocNeighbors_1.14.0  
#[101] cowplot_1.1.3          bitops_1.0-7           irlba_2.3.5.1          httpuv_1.6.14          patchwork_1.2.0       
#[106] R6_2.5.1               promises_1.2.1         KernSmooth_2.23-22     gridExtra_2.3          vipor_0.4.5           
#[111] parallelly_1.37.1      codetools_0.2-19       gtools_3.9.4           MASS_7.3-60            withr_3.0.0           
#[116] sctransform_0.3.5      GenomeInfoDbData_1.2.8 parallel_4.2.2         beachmat_2.12.0        grid_4.2.2            
#[121] tidyr_1.3.1            Rtsne_0.16             spatstat.explore_3.2-3 ggforce_0.4.1          shiny_1.8.0           
#[126] ggbeeswarm_0.7.2  









