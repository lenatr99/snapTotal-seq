
import scvelo as scv

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

adata = scv.read("HEK293T_snapTotal.h5ad")

adata.uns['cell_cycle_colors'] = ['#FF7F0E','#2CA02C','#6699FF']
adata.uns['seurat_clusters_colors'] = ['#f8766d','#a3a500','#00bf7d','#00b0f6','#e76bf3']

scv.pp.filter_genes(adata, min_shared_counts=30, min_cells_u=5)
scv.pp.normalize_per_cell(adata, enforce=True)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)

scv.pp.moments(adata, n_pcs=2, n_neighbors=20)

scv.tl.recover_dynamics(adata, fit_basal_transcription=True)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='pca', color='cell_cycle', alpha=0.4)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', alpha=0.4)

# velocity confidence score
scv.tl.velocity_confidence(adata)
keys = 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', basis='pca', vmin=0.4, vmax=1, size=300)

# infer latent time & velocity pseudotime
# as cell cycle is a cycle, instead of a linear process, an pre-defined endpoint is provided based on the velocity diagram
# based on the velocity diagram above, we chose the cell at the far-left as the pre-defined endpoint
# for any linear process, such prior information is not needed
# see more information at https://github.com/theislab/scvelo/issues

adata.obs['Endpoints'] = 0
adata.obs['Endpoints']['HEK293T_B37'] = 1

scv.tl.latent_time(adata,  end_key='Endpoints', min_likelihood=0.25)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', basis='pca', size=300)
scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot', basis='pca', size=300)

df = adata.obs
df.to_csv('HEK293T_scvelo_res.txt', sep='\t', index=True, header=True)
