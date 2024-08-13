library(Seurat)
library(SeuratDisk)
library(ggplot2)

seu <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqMastitisMilkBlood/Seurat/20240423_JEW_IntegratedSeurat_GranulocytesOnly_PhyloOrder.h5seurat')

dat <- data.frame(table(seu$phyloorder, seu$sample_ID))
colnames(dat) <- c('cluster', 'sample_ID', 'Freq')

dat$cluster <- factor(dat$cluster, levels = c(
  c('c30', 'c46', 'c29', 'c1', 'c33', 'c23', 'c18',
    'c32', 'c28', 'c26', 'c8', 'c5', 'c3', 'c10', 'c19',
    'c9', 'c24', 'c7', 'c6', 'c39', 'c16', 'c2', 
    'c11', 'c52', 'c58', 'c4', 'c60', 'c48', 'c12', 'c25')))
dat$sample_ID <- factor(dat$sample_ID, levels = rev(c('1312blood', '1630blood', '1634blood',
                                                      '1312milk', '1630milk', '1634milk')))


ggplot(dat, aes(fill=sample_ID, y=Freq, x=cluster)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = rev(c('indianred1', 'indianred3', 'indianred4',
                                   'dodgerblue1', 'dodgerblue3', 'dodgerblue4'))) +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


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

dat <- data.frame(table(seu$nodes, seu$sample_ID))
colnames(dat) <- c('nodes', 'sample_ID', 'Freq')

dat$sample_ID <- factor(dat$sample_ID, levels = c('1312blood', '1630blood', '1634blood',
                                                      '1312milk', '1630milk', '1634milk'))
dat2 <- subset(dat, nodes == 'node1')
ggplot(dat2, aes(x=1, y=Freq, fill=sample_ID))+ # pie chart for IPP cells
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_manual(values=c('indianred1', 'indianred3', 'indianred4',
                             'dodgerblue1', 'dodgerblue3', 'dodgerblue4')) +
  ggtitle('node1')

dat2 <- subset(dat, nodes == 'node2')
ggplot(dat2, aes(x=1, y=Freq, fill=sample_ID))+ # pie chart for IPP cells
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_manual(values=c('indianred1', 'indianred3', 'indianred4',
                             'dodgerblue1', 'dodgerblue3', 'dodgerblue4')) +
  ggtitle('node2')

dat2 <- subset(dat, nodes == 'node3')
ggplot(dat2, aes(x=1, y=Freq, fill=sample_ID))+ # pie chart for IPP cells
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_manual(values=c('indianred1', 'indianred3', 'indianred4',
                             'dodgerblue1', 'dodgerblue3', 'dodgerblue4')) +
  ggtitle('node3')
