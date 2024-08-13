library(Seurat)
library(future)
library(SeuratDisk)
library(ggplot2)

set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230129_JEW_IntegratedSeurat_AnnotatedClusters.h5seurat')
seu

Idents(seu) <- seu$cluster
DefaultAssay(seu) <- 'integrated'
seu <- BuildClusterTree(seu, 
                        dims = 1:30, 
                        assay = "PCA")
data.tree <- Tool(object = seu, 
                  slot = "BuildClusterTree") 
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, 
                               c('c53', 'c52', 'c3', 'c10', 'c5', 'c18', 'c1', 
                                 'c33', 'c19', 'c8', 'c26', 'c7', 'c29', 'c24', 
                                 'c9', 'c32', 'c6', 'c23', 'c16', 'c2', 'c39',
                                 'c46', 'c11', 'c28', 'c60', 'c12', 'c25', 'c48', 
                                 'c4', 'c58', 'c30', 'c59', 'c55', 'c51', 'c38', 
                                 'c43', 'c34', 'c20', 'c35', 'c57', 'c44', 'c41', 
                                 'c50', 'c21', 'c27', 'c36', 'c13', 'c31', 'c49',
                                 'c61', 'c17', 'c47', 'c40', 'c42', 'c22', 'c37', 
                                 'c14', 'c15', 'c45', 'c56', 'c54', 'c62'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)

Idents(seu) <- seu$cluster
levels(seu) <- c('c53', 'c52', 'c3', 'c10', 'c5', 'c18', 'c1', 
                'c33', 'c19', 'c8', 'c26', 'c7', 'c29', 'c24', 
                'c9', 'c32', 'c6', 'c23', 'c16', 'c2', 'c39',
                'c46', 'c11', 'c28', 'c60', 'c12', 'c25', 'c48', 
                'c4', 'c58', 'c30', 'c59', 'c55', 'c51', 'c38', 
                'c43', 'c34', 'c20', 'c35', 'c57', 'c44', 'c41', 
                'c50', 'c21', 'c27', 'c36', 'c13', 'c31', 'c49',
                'c61', 'c17', 'c47', 'c40', 'c42', 'c22', 'c37', 
                'c14', 'c15', 'c45', 'c56', 'c54', 'c62')
seu$phyloorder <- Idents(seu)

Idents(seu) <- seu$phyloorder
DefaultAssay(seu) <- 'RNA'
DotPlot(seu,
        features = c('PTPRC', # leukocyte
                     'GZMB', 'IL3RA', # pDC...also CD4+ CD14- MHCII-
                     'CSF3R', 'CXCL8', 'SRGN', # GRANULOCYTE
                     'LALBA', 'CSN2', # non-immune
                     'TYROBP', 'SIRPA', 'CD68', 'TLR4', 'ADGRE1', 'HCK', #pan-myeloid
                     'CD14', 'TREM2', 'AIF1', 'CSF1R', 'CST3', 'CD86', 'CD83', # monocyte/macrophage
                     'FCGR3A', # CD16 monocyte
                     'FLT3', 'IRF8', 'ENSBTAG00000009656', 'ENSBTAG00000013919', # cDC
                     'CD79A', 'CD79B', 'CD19', 'PAX5', 'MS4A1', 'BLNK', 'IRF4', 'JCHAIN', # B/ASC
                     'CD3E', 'CD3D', 'CD3G', 'CD247', 'ZAP70', 'CD2', 'ENSBTAG00000055197', 'CD4', 'CD8B', 'CD8A', # T
                     'NCR1', 'KLRB1', 'NKG2A', 'KLRD1', 'NKG7', 'GNLY', 'CTSW', 'PRF1', # ILC
                     'PCLAF', 'PCNA', 'UBE2C', 'MKI67', # cycling
                     'AHSP', 'HBB', 'HBM'), #RBC
        cols = c('khaki', 'firebrick'),
        group.by = 'phyloorder') & RotatedAxis()

ggsave('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Figures/DotPlot_AllCells_CanonicalGenes_FinalClusters_PhyloOrder.jpeg')

SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230130_JEW_IntegratedSeurat_AnnotatedClusters_PhyloOrder.h5seurat', overwrite = TRUE)
