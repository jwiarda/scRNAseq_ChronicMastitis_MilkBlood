library(Seurat)
library(future)
library(SeuratDisk)
library(clustree)
library(ggplot)
library(dplyr)
library(AUCell)
library(plyr)

set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly_PhyloOrder.h5seurat')
seu

# Identify gene signature ----
de <- read_xlsx("/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/FinalClustersDEGs_GranulocytesOnly.xlsx")
de <- subset(de, avg_log2FC > 0) 
Idents(seu) <- seu$phyloorder
de30 <- subset(de, cluster == 'c30')
de46 <- subset(de, cluster == 'c46')
node1DE <- Reduce(intersect, list(de30$gene, de46$gene))
node1DE <- sort(node1DE)
length(node1DE)

# Perform gene set enrichment analysis ----
geneSets <- list(node1DE = node1DE)
exprMatrix <- as.matrix(seu[["RNA"]]@counts) 
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) # calculate cell signature AUC score for each gene set in each cell

clusterID <- as.data.frame(seu[["phyloorder"]])
clusterID <- t(clusterID)
AUCs <- as.data.frame(getAUC(cells_AUC))
AUCs <- rbind(AUCs, clusterID)
AUCs <- t(AUCs)
AUCs <- as.data.frame(AUCs)
head(AUCs)

geneSetName <- rownames(cells_AUC)[grep("node1DE", rownames(cells_AUC))]
thres = 0
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=thres)

coords <- Embeddings(seu[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$node1DE)
coords$`AUCs$node1DE` <- as.numeric(coords$`AUCs$node1DE`)
ggplot() + # make a tsne plot
  geom_point(data=coords, aes(x=tSNE_1, y=tSNE_2, color=coords[,3])) +
  scale_color_gradientn(colors = c('dodgerblue4', 'slateblue4', 'violetred', 'orange', 'gold')) +
  theme_void()

AUCs$node1DE <- as.numeric(AUCs$node1DE)
seu$node1Sig <- AUCs$node1DE
VlnPlot(seu, features = 'node1Sig', pt.size = 0.001)+
  stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
  scale_fill_manual(values=rep(c('grey80'), times = 30))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df <- data.frame(seu$phyloorder, seu$node1Sig)
colnames(df) <- c('cluster', 'AUC')
dat <- data_summary(df, varname="AUC", 
                    groupnames=c("cluster"))

dat <- dat[order(dat$AUC), ]
clusOrder <- dat$cluster
seu$cluster <- factor(seu$cluster, 
                      levels=rev(c(clusOrder)))

VlnPlot(seu, features = 'node1Sig', pt.size = 0.001, group.by = 'cluster')+
  stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
  scale_fill_manual(values=rep(c('grey80'), times = 30))

# Stacked bar of tissue type proportions ----
prop <- data.frame(table(seu$phyloorder, seu$tissue))
colnames(prop) <- c('cluster', 'tissue', 'Freq')

prop$cluster <- factor(prop$cluster, levels =c(clusOrder))
prop$tissue <- factor(prop$tissue, levels = rev(c('blood', 'milk')))


ggplot(prop, aes(fill=tissue, y=Freq, x=cluster)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = rev(c('indianred3', 'dodgerblue3'))) +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


dat2 <- data.frame(prop.table(table(seu$cluster, seu$tissue), margin = 1))
dat2 <- subset(dat2, Var2 == 'milk')
dat2 <- dat2[order(dat2$Freq), ]
clusOrder2 <- dat2$Var1

dat$cluster <- factor(dat$cluster, levels =c(clusOrder2))
ggplot(dat, aes(x = cluster, y = 1, fill = AUC))+ 
  geom_tile(color = 'black', size = .2)+
  scale_fill_gradientn(colors = c('gold', 'orange', 'red', 'red3'), limits = c(min(dat$AUC), max(dat$AUC))) + 
  theme_classic()

FeaturePlot(seu, reduction = 'tsne', 
            features = 'node1Sig', 
            cols = c('gold', 'orange', 'red', 'red3'), 
            min.cutoff = (min(dat$AUC)),
            max.cutoff = max(dat$AUC))

# Scatterplot ----
seu$cluster <- factor(seu$cluster, 
                      levels=rev(c(clusOrder2)))
VlnPlot(seu, 
        group.by = 'cluster', 
        split.by = 'tissue', 
        features = 'node1Sig', 
        cols = c('indianred3', 'dodgerblue3'),
        #split.plot = TRUE,
        pt.size = 0.001)

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.4 LTS

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
#  [1] plyr_1.8.9            AUCell_1.18.1         clustree_0.5.0        ggraph_2.1.0          stringr_1.5.1         dplyr_1.1.4          
#[7] scales_1.3.0          ggplot2_3.5.0         readxl_1.4.3          writexl_1.4.2         SeuratDisk_0.0.0.9020 future_1.33.1        
#[13] SeuratObject_4.1.3    Seurat_4.3.0.1       

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.4                  R.utils_2.12.2              spatstat.explore_3.2-3      reticulate_1.35.0          
#[5] tidyselect_1.2.0            RSQLite_2.3.1               AnnotationDbi_1.58.0        htmlwidgets_1.6.4          
#[9] grid_4.2.2                  BiocParallel_1.30.4         Rtsne_0.16                  munsell_0.5.0              
#[13] ScaledMatrix_1.4.1          codetools_0.2-19            ica_1.0-3                   miniUI_0.1.1.1             
#[17] withr_3.0.0                 spatstat.random_3.1-6       colorspace_2.1-0            progressr_0.14.0           
#[21] Biobase_2.56.0              rstudioapi_0.15.0           stats4_4.2.2                SingleCellExperiment_1.18.1
#[25] ROCR_1.0-11                 tensor_1.5                  listenv_0.9.1               labeling_0.4.3             
#[29] MatrixGenerics_1.8.1        GenomeInfoDbData_1.2.8      polyclip_1.10-4             bit64_4.0.5                
#[33] farver_2.1.1                parallelly_1.37.1           vctrs_0.6.5                 generics_0.1.3             
#[37] R6_2.5.1                    GenomeInfoDb_1.32.4         ggbeeswarm_0.7.2            graphlayouts_1.0.0         
#[41] rsvd_1.0.5                  locfit_1.5-9.8              miloR_1.4.0                 hdf5r_1.3.8                
#[45] bitops_1.0-7                spatstat.utils_3.0-3        cachem_1.0.8                DelayedArray_0.22.0        
#[49] promises_1.2.1              beeswarm_0.4.0              gtable_0.3.4                beachmat_2.12.0            
#[53] globals_0.16.2              goftest_1.2-3               tidygraph_1.2.3             rlang_1.1.3                
#[57] splines_4.2.2               lazyeval_0.2.2              spatstat.geom_3.2-5         BiocManager_1.30.22        
#[61] reshape2_1.4.4              abind_1.4-5                 httpuv_1.6.14               tools_4.2.2                
#[65] ellipsis_0.3.2              RColorBrewer_1.1-3          BiocGenerics_0.42.0         ggridges_0.5.4             
#[69] Rcpp_1.0.12                 sparseMatrixStats_1.8.0     zlibbioc_1.42.0             purrr_1.0.2                
#[73] RCurl_1.98-1.12             deldir_1.0-9                pbapply_1.7-2               viridis_0.6.4              
#[77] cowplot_1.1.3               S4Vectors_0.34.0            zoo_1.8-12                  SummarizedExperiment_1.26.1
#[81] ggrepel_0.9.5               cluster_2.1.4               magrittr_2.0.3              data.table_1.15.2          
#[85] scattermore_1.2             lmtest_0.9-40               RANN_2.6.1                  fitdistrplus_1.1-11        
#[89] matrixStats_1.0.0           patchwork_1.2.0             mime_0.12                   xtable_1.8-4               
#[93] XML_3.99-0.14               IRanges_2.30.1              gridExtra_2.3               compiler_4.2.2             
#[97] tibble_3.2.1                KernSmooth_2.23-22          crayon_1.5.2                R.oo_1.25.0                
#[101] htmltools_0.5.7             later_1.3.2                 tidyr_1.3.1                 DBI_1.2.2                  
#[105] tweenr_2.0.2                MASS_7.3-60                 Matrix_1.6-1                cli_3.6.2                  
#[109] R.methodsS3_1.8.2           parallel_4.2.2              igraph_1.5.1                GenomicRanges_1.48.0       
#[113] pkgconfig_2.0.3             sp_2.0-0                    plotly_4.10.4               spatstat.sparse_3.0-2      
#[117] annotate_1.74.0             vipor_0.4.5                 XVector_0.36.0              digest_0.6.34              
#[121] sctransform_0.3.5           RcppAnnoy_0.0.21            graph_1.74.0                spatstat.data_3.0-1        
#[125] Biostrings_2.64.1           cellranger_1.1.0            leiden_0.4.3                uwot_0.1.16                
#[129] edgeR_3.38.4                DelayedMatrixStats_1.18.2   GSEABase_1.58.0             shiny_1.8.0                
#[133] gtools_3.9.4                lifecycle_1.0.4             nlme_3.1-163                jsonlite_1.8.8             
#[137] BiocNeighbors_1.14.0        viridisLite_0.4.2           limma_3.52.4                fansi_1.0.6                
#[141] pillar_1.9.0                lattice_0.21-8              ggrastr_1.0.2               KEGGREST_1.36.3            
#[145] fastmap_1.1.1               httr_1.4.7                  survival_3.5-7              glue_1.7.0                 
#[149] png_0.1-8                   bit_4.0.5                   ggforce_0.4.1               stringi_1.8.3              
#[153] blob_1.2.4                  BiocSingular_1.12.0         memoise_2.0.1               irlba_2.3.5.1              
#[157] future.apply_1.11.1 