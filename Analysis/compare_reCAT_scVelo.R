# Compare the cell cycle trajectory inferred by reCAT and scVelo
# As cell cycle is a cycle, reCAT and scVelo may select different cells
# as the start point (and end point). Therefore, we use the best match between 
# two trajectories to determine the consistency between two methods.

library(ggplot2)
rm(list = ls())

# load scVelo and reCAT results
velo_res = read.table('HEK293T_scvelo_res.txt', header = T, 
                      row.names = 1, sep = '\t', as.is = T)
# order cells by latent time
velo_res = velo_res[order(velo_res$latent_time), ]

reCAT_df = read.table('293T_reCAT_result.txt', 
                      header = T, sep = '\t', as.is = T)
rownames(reCAT_df) = reCAT_df$cell_id
# cells have been ordered along cell cycle in reCAT input file

ncell = nrow(reCAT_df)

corr = c()
# change the start point and calculate correlation coefficients 
# to find the best match
for(i in 2:ncell){
  tmp_df = reCAT_df[c(i:ncell, 1:(i-1)), ]
  tmp_df$ord = 1:ncell
  corr = c(corr, cor(1:ncell, tmp_df[rownames(velo_res), 'ord']))
}
reCAT_df$ord = 1:ncell
corr = c(cor(1:ncell, reCAT_df[rownames(velo_res), 'ord']), corr)

# best match
print(max(corr))
best_match = which(corr == max(corr))

# reorder with the start point giving the best match
reCAT_df = reCAT_df[c(best_match:ncell, 1:(best_match-1)), ]
reCAT_df$ord = 1:ncell

# plot
plot_df = data.frame('reCAT' = reCAT_df[rownames(velo_res), 'ord'],
                     'scVelo' = 1:ncell)

ggplot(data = plot_df, mapping = aes(x = reCAT, y = scVelo)) + 
  geom_point(color = 'gray60') + 
  geom_smooth(method='lm') + 
  theme_classic() + 
  labs(title="",x="Order by reCAT", y = "Order by latent time") + 
  theme(axis.text.x = element_text(colour="black",size=12,face = "plain"),
        axis.text.y = element_text(colour = "black", size = 12, face = "plain"),
        axis.title.x = element_text(colour = "black", size = 16),
        axis.title.y = element_text(colour = "black", size = 16),
        aspect.ratio = 1)

