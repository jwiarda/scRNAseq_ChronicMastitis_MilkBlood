library(Seurat)
library(future)
library(SeuratDisk)
library(ggplot2)
library(writexl)
library(readxl)
library(scales)
library(dplyr)
library(stringr) 

set.seed(5)

# setup plan for mulitprocessing
if (future::supportsMulticore()){
  future::plan(multicore, workers=4)
} else {
  future::plan(multisession, workers=4)
}

# 200GB limit
options(future.globals.maxSize = 200000 * 1024^2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230130_JEW_IntegratedSeurat_AnnotatedClusters_PhyloOrder.h5seurat')
seu

# Cluster vs all ----
DefaultAssay(seu) <- 'RNA'
Idents(seu) <- seu$phyloorder

de <- FindAllMarkers(seu)
de <- subset(de, p_val_adj < 0.05) # make sure the adjusted p-values are still < 0.05 since some genes in DE list have p_val_adj > 0.05

# Save DE gene list:
write_xlsx(x = de, 
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/FinalClustersDEGs_AllCells.xlsx",
           col_names = TRUE)

# Save background gene list:
bg <- data.frame(rownames(seu))
colnames(bg) <- 'gene'
write_xlsx(x = bg, 
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/BackgroundGeneList_AllCells.xlsx",
           col_names = TRUE)

# Make heatmap of top DEGs:
#de <- read_excel('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/FinalClustersDEGs.xlsx')
de <- subset(de, avg_log2FC > 0) # only take genes enriched in the clusters 
Idents(seu) <- seu$phyloorder
topgenes <- de %>% group_by(cluster) %>% top_n(5, avg_log2FC) # only plot top 5 genes per cluster, as determined by highest average logFC values
DoHeatmap(subset(seu, downsample = 100), # take only 50 cells per cluster for plotting
          features = as.character(topgenes$gene), 
          assay = "RNA", 
          group.colors = c('deepskyblue', rep('chartreuse4', 30), 'grey80', 
                           rep('violetred', 8), 'goldenrod1', rep('slateblue', 7),
                           rep('chocolate2', 12), 'slateblue', 'violetred'),
          disp.min = -1, 
          disp.max = 2) +
  scale_fill_gradientn(colors = c("darkturquoise", "grey90", "indianred1", "red"))

# Perform pairwise DGE analysis ----

Idents(seu) <- seu$phyloorder
clusters <- unique(Idents(seu)) # identify all of our cluster IDs
pairwise <- combn(clusters, 2) # create all pairwise combinations of cluster IDs
p1 <- pairwise[1,] 
p2 <- pairwise[2,] 
comps1 <- data.frame(p1, p2)
colnames(comps1) <- c('pop1', 'pop2')
comps2 <- data.frame(p2, p1)
colnames(comps2) <- c('pop1', 'pop2')
allcomps <- rbind(comps1, comps2)
dim(allcomps) # see how many comparisons are to be made, then break them up so we don't have any session time out issues

comps <- allcomps[1:500,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_1to500.xlsx",
           col_names = TRUE)

comps <- allcomps[501:1000,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_501to1000.xlsx",
           col_names = TRUE)


comps <- allcomps[1001:1500,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_1001to1500.xlsx",
           col_names = TRUE)

comps <- allcomps[1501:2000,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_1501to2000.xlsx",
           col_names = TRUE)

comps <- allcomps[2001:2500,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_2001to2500.xlsx",
           col_names = TRUE)

comps <- allcomps[2501:3000,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_2501to3000.xlsx",
           col_names = TRUE)

comps <- allcomps[3001:3500,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_3001to3500.xlsx",
           col_names = TRUE)

comps <- allcomps[3501:3782,]
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
           path = "/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_3501to3782.xlsx",
           col_names = TRUE)

# Merge all de into one file:
de1 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_1to500.xlsx')
de2 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_501to1000.xlsx')
de3 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_1001to1500.xlsx')
de4 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_1501to2000.xlsx')
de5 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_2001to2500.xlsx')
de6 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_2501to3000.xlsx')
de7 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_3001to3500.xlsx')
de8 <- read_xlsx('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_pairwise_3501to3782.xlsx')
de <- bind_rows(de1, de2, de3, de4, de5, de6, de7, de8)
write.table(de, file='/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/DGE/AllCellsDEGs_AllPairwiseDEGs.tsv', quote=FALSE, sep='\t', col.names = NA)

# Find total number of DEGs for each comparison
sum <- data.frame(table(de$comparison))
colnames(sum) <- c('comparison', 'count')
sum$pop1 <- sub("\\v.*", "", sum$comparison)
sum$pop2 <- sub('.*v', '', sum$comparison)

sum$pop1 <- factor(sum$pop1,levels = c(levels(seu$phyloorder)))
sum$pop2 <- factor(sum$pop2,levels = c(levels(seu$phyloorder)))

ggplot(sum, aes(pop1, pop2, fill = count)) +
  geom_tile()+
  #scale_fill_gradientn(colours = c('slateblue4', 'violetred', 'orange', 'gold'))+ 
  theme_classic()+
  scale_fill_gradientn(colours = c('beige', 'gold', 'orange','red', 'darkred', 'brown', 'black'),
                       limits = c(0, 4100), oob=squish) +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
