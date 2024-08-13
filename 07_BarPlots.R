library(Seurat)
library(SeuratDisk)
library(ggplot2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20230130_JEW_IntegratedSeurat_AnnotatedClusters_PhyloOrder.h5seurat')

dat <- data.frame(table(seu$celltype, seu$sample_ID))
colnames(dat) <- c('celltype', 'sample_ID', 'Freq')

dat$celltype <- factor(dat$celltype, levels = rev(c('granulocyte', 'monocyte/macrophage/cDC', 
                                                     'B/ASC', 'T/ILC', 'pDC', 'non-immune')))
dat$sample_ID <- factor(dat$sample_ID, levels = c('1312blood', '1630blood', '1634blood',
                                                  '1312milk', '1630milk', '1634milk'))
                    
ggplot(dat, aes(fill=celltype, y=Freq, x=sample_ID)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = rev(c('chartreuse4', 'violetred', 'slateblue', 'chocolate2', 
                               'deepskyblue2', 'goldenrod1'))) +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dat <- data.frame(table(seu$phyloorder, seu$sample_ID))
colnames(dat) <- c('cluster', 'sample_ID', 'Freq')

dat$cluster <- factor(dat$cluster, levels = c(
  c('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 
    'c12', 'c16', 'c18', 'c19', 'c23', 'c24', 'c25', 'c26', 
    'c28', 'c29', 'c30', 'c32', 'c33', 'c39', 'c46', 'c48', 'c52', 'c58', 'c60',
    
    'c20', 'c34', 'c35', 'c38', 'c43', 'c51', 'c55', 'c57', 'c59', 'c62',
    
    'c13', 'c21', 'c27', 'c31', 'c36', 'c41', 'c50', 'c54',
    
    'c14', 'c15', 'c17', 'c22', 'c37', 'c40', 'c42', 'c45', 'c47', 'c49', 'c56', 'c61',
    
    'c53',
    
    'c44')))
dat$sample_ID <- factor(dat$sample_ID, levels = rev(c('1312blood', '1630blood', '1634blood',
                                                  '1312milk', '1630milk', '1634milk')))


ggplot(dat, aes(fill=sample_ID, y=Freq, x=cluster)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = rev(c('indianred1', 'indianred3', 'indianred4',
                                   'dodgerblue1', 'dodgerblue3', 'dodgerblue4'))) +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dat <- data.frame(table(seu$celltype, seu$sample_ID))
colnames(dat) <- c('celltype', 'sample_ID', 'Freq')
dat2 <- data.frame(table(seu$sample_ID))
colnames(dat2) <- c('sample_ID', 'Freq')
dat2$celltype <- rep('TotalCells', ncol(dat2))
dat2 <- dat2[, c('celltype', 'sample_ID', 'Freq')] # leave the row index blank to keep all rows
dat <- rbind(dat, dat2)


dat$celltype <- factor(dat$celltype, levels = rev(c('granulocyte', 'monocyte/macrophage/cDC', 
                                                    'B/ASC', 'T/ILC', 'pDC', 'non-immune', 'TotalCells')))
dat$sample_ID <- factor(dat$sample_ID, levels = rev(c('1312blood', '1630blood', '1634blood',
                                                      '1312milk', '1630milk', '1634milk')))


ggplot(dat, aes(fill=sample_ID, y=Freq, x=celltype)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = rev(c('indianred1', 'indianred3', 'indianred4',
                                   'dodgerblue1', 'dodgerblue3', 'dodgerblue4'))) +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

