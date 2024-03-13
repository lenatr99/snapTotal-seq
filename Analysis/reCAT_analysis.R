# reCAT analysis
# reCAT package can be downloaded from https://github.com/tinglab/reCAT
rm(list=ls())

# Optional: set working directory
#setwd('[WORKING DIRECTORY]')

# using exon data for reCAT analysis
exon_df = read.table('GSE202126_HEK293T_exon_UMI_count_matrix.txt',
                     header = T, sep = '\t', row.names = 1, as.is = T)
nucidx = grep('^MT-', exon_df$gene_symbol, invert = T)
exon_df = exon_df[nucidx,]
exon_df = exon_df[!duplicated(exon_df$gene_symbol),]
n = dim(exon_df)[2]
exon_mat = as.matrix(exon_df[, 2:n])
dim(exon_mat)

# remove lowly exp genes
whether_exp = exon_mat > 1
keep = (rowSums(whether_exp) >= 5)
table(keep)
exon_mat = exon_mat[keep, ]
gene_list = rownames(exon_mat)
ngene = length(gene_list)
cell_list = colnames(exon_mat)
ncell = length(cell_list)
rownames(exon_mat) = exon_df[gene_list, 'gene_symbol']
exon_total_umi = colSums(exon_mat)
avg_exon_umi = median(exon_total_umi)

# normalization
exon_norm_mat = matrix(0, nrow = ngene, ncol = ncell)
rownames(exon_norm_mat) = rownames(exon_mat)
colnames(exon_norm_mat) = colnames(exon_mat)
for(i in 1:ncell)
{
  exon_norm_mat[,i] = exon_mat[,i] / exon_total_umi[i] * avg_exon_umi
}
exon_norm_mat = log2(exon_norm_mat + 1)

# run reCAT
# Please refer to https://github.com/tinglab/reCAT for more detailed instructions. 

# Optional: change to the directory with reCAT scripts to load reCAT functions
source("get_test_exp.R")
test_exp <- get_test_exp(exon_norm_mat)

source("get_ordIndex.R")
ordIndex <- get_ordIndex(test_exp, 1)

source("get_score.R")
score_result <- get_score(t(test_exp))

source("plot.R")
plot_bayes(score_result$bayes_score, ordIndex)
plot_mean(score_result$mean_score, ordIndex)

source("get_hmm.R")
myord = c(90:1, 143:91)
rd = data.frame('start' = c(71,21,110),
                'end' = c(82,40,125))
rd = as.matrix(rd)
hmm_result <- get_hmm_order(bayes_score = score_result$bayes_score, 
                            mean_score = score_result$mean_score, 
                            ordIndex = ordIndex, cls_num = 3, 
                            myord = myord, rdata = rd)
hmm_result

hmm_df = data.frame('cell_id' = rownames(score_result$mean_score)[ordIndex][myord],
                    'cell_cycle' = hmm_result)

# output result to working directory
write.table(hmm_df, '[WORKING DIRECTORY]/293T_reCAT_result.txt', col.names = T, row.names = F, 
            sep = '\t', quote = F)



