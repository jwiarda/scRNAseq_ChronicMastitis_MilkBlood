library(Seurat)
library(future)
library(SeuratDisk)
library(ggplot2)
library(writexl)
library(readxl)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly_PhyloOrder.h5seurat')
Idents(seu) <- seu$phyloorder
n1 <- rep('node1', 2)
n2 <- rep('node2', 22)
n3 <- rep('node3', 6)
nodes <- c(n1, n2, n3)
seu$nodes <- seu$phyloorder
Idents(seu) <- seu$nodes
names(nodes) <- levels(seu)
seu <- RenameIdents(seu, nodes)
seu$nodes <- Idents(seu)

#########################################
seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly_PhyloOrder.h5seurat')
DefaultAssay(seu) <- 'RNA'

de <- read_xlsx("/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/FinalClustersDEGs_GranulocytesOnly.xlsx")
de <- subset(de, avg_log2FC > 0) 
Idents(seu) <- seu$phyloorder

de30 <- subset(de, cluster == 'c30')
de46 <- subset(de, cluster == 'c46')
node1DE <- Reduce(intersect, list(de30$gene, de46$gene))
node1DE <- sort(node1DE)
length(node1DE)

DotPlot(seu, features = c(node1DE), cols = c('goldenrod1', 'firebrick')) + RotatedAxis() # this dot plot has so many genes, let's just show a few with highest logFC...
DotPlot(seu, features = c(node1DE[1:105]), cols = c('goldenrod1', 'firebrick'), col.min = -1.5, col.max = 2) + RotatedAxis() # this dot plot has so many genes, let's just show a few with highest logFC...
DotPlot(seu, features = c(node1DE[106:210]), cols = c('goldenrod1', 'firebrick'), col.min = -1.5, col.max = 2) + RotatedAxis() # this dot plot has so many genes, let's just show a few with highest logFC...


de30 <- subset(de, cluster == 'c30' & avg_log2FC > 1)
de46 <- subset(de, cluster == 'c46' & avg_log2FC > 1)
node1DE <- Reduce(intersect, list(de30$gene, de46$gene))
node1DE <- sort(node1DE)
length(node1DE)

DotPlot(seu, features = c(node1DE), cols = c('goldenrod1', 'firebrick'), col.min = -1, col.max = 2) + RotatedAxis() # this dot plot has so many genes, let's just show a few with highest logFC...


de19 <- subset(de, cluster == 'c19')
de60 <- subset(de, cluster == 'c60')
de30 <- subset(de, cluster == 'c30')
de46 <- subset(de, cluster == 'c46')
node1DE <- Reduce(intersect, list(de30$gene, de46$gene, de19$gene, de60$gene))
node1DE <- sort(node1DE)
length(node1DE)

dot <- DotPlot(seu, features = c(node1DE), cols = c('goldenrod1', 'firebrick'), col.min = -1, col.max = 2) + RotatedAxis() # this dot plot has so many genes, let's just show a few with highest logFC...
dot
dot$data # extract dot plot data metrics
tab <- data.frame(dot$data)

dat <- data.frame(prop.table(table(seu$cluster, seu$tissue), margin = 1))
dat <- subset(dat, Var2 == 'milk')
dat <- dat[order(dat$Freq), ]
clusOrder <- dat$Var1

dat <- data.frame(table(seu$phyloorder, seu$tissue))
colnames(dat) <- c('cluster', 'tissue', 'Freq')

dat$cluster <- factor(dat$cluster, levels =c(clusOrder))
dat$tissue <- factor(dat$tissue, levels = rev(c('blood', 'milk')))


ggplot(dat, aes(fill=tissue, y=Freq, x=cluster)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = rev(c('indianred3', 'dodgerblue3'))) +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Idents(seu) <- seu$cluster
seu$cluster <- factor(seu$cluster, 
                            levels=rev(c(clusOrder)))
DotPlot(seu, features = c(node1DE), cols = c('goldenrod1', 'firebrick'), col.min = -1, col.max = 2) + RotatedAxis() # this dot plot has so many genes, let's just show a few with highest logFC...

sessionInfo()
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
#[4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggplot2_3.5.0         readxl_1.4.3          writexl_1.4.2         SeuratDisk_0.0.0.9020
#[5] future_1.33.1         SeuratObject_4.1.3    Seurat_4.3.0.1       

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-9           ellipsis_0.3.2        
#[5] ggridges_0.5.4         rstudioapi_0.15.0      spatstat.data_3.0-1    leiden_0.4.3          
#[9] listenv_0.9.1          farver_2.1.1           ggrepel_0.9.5          bit64_4.0.5           
#[13] fansi_1.0.6            codetools_0.2-19       splines_4.2.2          polyclip_1.10-4       
#[17] jsonlite_1.8.8         ica_1.0-3              cluster_2.1.4          png_0.1-8             
#[21] uwot_0.1.16            shiny_1.8.0            sctransform_0.3.5      spatstat.sparse_3.0-2 
#[25] compiler_4.2.2         httr_1.4.7             Matrix_1.6-1           fastmap_1.1.1         
#[29] lazyeval_0.2.2         cli_3.6.2              later_1.3.2            htmltools_0.5.7       
#[33] tools_4.2.2            igraph_1.5.1           gtable_0.3.4           glue_1.7.0            
#[37] RANN_2.6.1             reshape2_1.4.4         dplyr_1.1.4            Rcpp_1.0.12           
#[41] scattermore_1.2        cellranger_1.1.0       vctrs_0.6.5            spatstat.explore_3.2-3
#[45] nlme_3.1-163           progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.1-6 
#[49] stringr_1.5.1          globals_0.16.2         mime_0.12              miniUI_0.1.1.1        
#[53] lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3          MASS_7.3-60           
#[57] zoo_1.8-12             scales_1.3.0           promises_1.2.1         spatstat.utils_3.0-3  
#[61] parallel_4.2.2         RColorBrewer_1.1-3     reticulate_1.35.0      pbapply_1.7-2         
#[65] gridExtra_2.3          stringi_1.8.3          rlang_1.1.3            pkgconfig_2.0.3       
#[69] matrixStats_1.0.0      lattice_0.21-8         ROCR_1.0-11            purrr_1.0.2           
#[73] tensor_1.5             patchwork_1.2.0        htmlwidgets_1.6.4      labeling_0.4.3        
#[77] cowplot_1.1.3          bit_4.0.5              tidyselect_1.2.0       parallelly_1.37.1     
#[81] RcppAnnoy_0.0.21       plyr_1.8.9             magrittr_2.0.3         R6_2.5.1              
#[85] generics_0.1.3         pillar_1.9.0           withr_3.0.0            fitdistrplus_1.1-11   
#[89] survival_3.5-7         abind_1.4-5            sp_2.0-0               tibble_3.2.1          
#[93] future.apply_1.11.1    crayon_1.5.2           hdf5r_1.3.8            KernSmooth_2.23-22    
#[97] utf8_1.2.4             spatstat.geom_3.2-5    plotly_4.10.4          grid_4.2.2            
#[101] data.table_1.15.2      digest_0.6.34          xtable_1.8-4           tidyr_1.3.1           
#[105] httpuv_1.6.14          munsell_0.5.0          viridisLite_0.4.2 