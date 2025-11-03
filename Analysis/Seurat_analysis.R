library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(patchwork)
library(ggplot2)
rm(list=ls())

# optional: set the working directory
setwd('/Users/lenatrnovec/scTotalRNA/snapTotal-seq/Analysis')

# load data
exon_df = read.table('DICTY_exon_UMI_count_matrix.txt',
                     header = T, sep = '\t', row.names = 2, as.is = T)
nucidx = grep('^MT-', exon_df$gene_id, invert = T)
exon_df = exon_df[nucidx,]
exon_df = exon_df[!duplicated(exon_df$gene_id),]
n = dim(exon_df)[2]
exon_mat = as.matrix(exon_df[, 2:n])
dim(exon_mat)

# remove lowly exp genes
whether_exp = exon_mat > 0
keep = (rowSums(whether_exp) >= 0)
table(keep)
exon_mat = exon_mat[keep, ]
gene_list = rownames(exon_mat)
cell_list = colnames(exon_mat)
rownames(exon_mat) = exon_df[gene_list, 'gene_id']
exon_total_umi = colSums(exon_mat)

# Optional: load cell cycle assignment from reCAT
# reCAT_df = read.table('293T_reCAT_result.txt', header = T,
#                       sep = '\t', row.names = 1, as.is = T)
# head(reCAT_df)


# run Seurat
d = CreateSeuratObject(counts = exon_mat, project = "293T",
                       min.cells = 5, min.features = 2000)
d = NormalizeData(d, normalization.method = "LogNormalize", 
                  scale.factor = median(exon_total_umi))
d = FindVariableFeatures(d, selection.method = "vst", nfeatures = 500)
VariableFeaturePlot(d)




all.genes = rownames(d)
d = ScaleData(d, features = all.genes)

# Optional: add cell cycle assignment to Seurat object
# d@meta.data[['cell_cycle']] = reCAT_df[cell_list, 1]

vf <- VariableFeatures(d)

# pick a safe number of PCs given tiny sample size
safe_npcs <- max(2, min(ncol(d), length(vf)) - 1)  # here: 5 cells -> 4 PCs
safe_npcs
# [1] 4

d <- RunPCA(d, features = vf, npcs = safe_npcs)
DimPlot(d, reduction = 'pca')
DimPlot(d, reduction = 'pca', group.by = 'cell_cycle') + 
  scale_color_manual(values = c('#FF7F0E','#2CA02C','#6699FF'))

# UMAP and unsupervised clustering
ElbowPlot(d)
d = RunUMAP(d, dims = 1:6)
d = FindNeighbors(d, dims = 1:6, k.param = 10)
d = FindClusters(d, resolution = 0.4)
DimPlot(d, reduction = 'umap')

# generate input files for scVelo
# Load intron data
intron_df = read.table('DICTY_intron_UMI_count_matrix.txt',
                       header=T, row.names = 2, sep='\t', as.is = T)
intron_mat = as.matrix(intron_df[gene_list, cell_list])
dim(intron_mat)
rownames(intron_mat) = exon_df[gene_list, 'gene_id']

# Add unspliced data to Seurat object
d[["unspliced"]] = CreateAssayObject(counts = intron_mat)
d[["spliced"]] = CreateAssayObject(counts = exon_mat)
d[["RNA"]] = CreateAssayObject(counts = exon_mat)

# Output
SaveH5Seurat(d, filename = "DICTY_snapTotal.h5Seurat", overwrite = T)
Convert("DICTY_snapTotal.h5Seurat", dest = "h5ad", overwrite = T)



