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

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly.h5seurat')
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
                               rev(c('c30', 'c46', 'c29', 'c1', 'c33', 'c23', 'c18',
                                 'c32', 'c28', 'c26', 'c8', 'c5', 'c3', 'c10', 'c19',
                                 'c9', 'c24', 'c7', 'c6', 'c39', 'c16', 'c2', 
                                 'c11', 'c52', 'c58', 'c4', 'c60', 'c48', 'c12', 'c25')))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)

Idents(seu) <- seu$cluster
levels(seu) <- c('c30', 'c46', 'c29', 'c1', 'c33', 'c23', 'c18',
                 'c32', 'c28', 'c26', 'c8', 'c5', 'c3', 'c10', 'c19',
                 'c9', 'c24', 'c7', 'c6', 'c39', 'c16', 'c2', 
                 'c11', 'c52', 'c58', 'c4', 'c60', 'c48', 'c12', 'c25')
seu$phyloorder <- Idents(seu)

Idents(seu) <- seu$phyloorder
DefaultAssay(seu) <- 'RNA'

SaveH5Seurat(seu, '/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly_PhyloOrder.h5seurat', overwrite = TRUE)
