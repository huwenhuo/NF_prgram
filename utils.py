import scanpy as sc

def preprocess_and_transfer(adata, n_top_genes=2000, n_neighbors=15, n_pcs=30, leiden_res=0.6):
    """
    Preprocess an AnnData object using HVGs, run PCA, neighbors, UMAP, and Leiden clustering.
    Transfers the embeddings and clustering back to the original object.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object (original).
    n_top_genes : int
        Number of highly variable genes to select.
    n_neighbors : int
        Number of neighbors for neighbor graph.
    n_pcs : int
        Number of principal components for PCA and neighbors.
    leiden_res : float
        Resolution for Leiden clustering.
        
    Returns
    -------
    adata_hvg : AnnData
        HVG-subset AnnData with PCA, neighbors, UMAP, and Leiden.
    adata : AnnData
        Original AnnData with transferred PCA, UMAP, neighbors, and Leiden.
    """
    # Normalize and log-transform the original object
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat_v3')
    
    # Subset to HVGs and make a copy
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    
    # PCA
    sc.tl.pca(adata_hvg, svd_solver='arpack')
    
    # Neighbors
    sc.pp.neighbors(adata_hvg, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # UMAP
    sc.tl.umap(adata_hvg)
    
    # Leiden clustering
    sc.tl.leiden(adata_hvg, resolution=leiden_res)
    
    # Transfer results back to original object
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
    adata.uns['neighbors'] = adata_hvg.uns['neighbors']
    adata.obs['leiden'] = adata_hvg.obs['leiden']
    
    return adata_hvg, adata

import scanpy as sc
import pandas as pd

def prepare_dotplot_data(adata, genes, groupby='celltype_sample'):
    """
    Prepare dotplot data (mean expression and fraction expressing) for custom plotting.

    Parameters
    ----------
    adata : AnnData
        AnnData object.
    genes : list of str
        List of gene names to include.
    groupby : str
        Column in adata.obs to group by (e.g., cell_type, sample, or combined).

    Returns
    -------
    df_plot : pd.DataFrame
        Tidy DataFrame with columns: group, gene, avg_expr, frac_expr
    """
    # Ensure genes is a list
    if isinstance(genes, str):
        genes = [genes]

    # Aggregate counts
    aggregated = sc.get.aggregate(adata[:, genes], by=groupby, func=["sum", "count_nonzero"])

    # Convert to DataFrames and reset index
    sum_df = aggregated.to_df(layer="sum").reset_index().rename(columns={"index": groupby})
    count_df = aggregated.to_df(layer="count_nonzero").reset_index().rename(columns={"index": groupby})

    # Compute mean in expressing cells
    mean_df = sum_df.copy()
    mean_df[genes] = sum_df[genes] / count_df[genes]

    # Compute fraction expressing: number of expressing cells / total cells in group
    total_cells_per_group = adata.obs[groupby].value_counts()
    frac_df = count_df.copy()
    for gene in genes:
        frac_df[gene] = count_df[gene] / count_df[groupby].map(total_cells_per_group)

    # Melt to tidy format
    df_mean = mean_df.melt(id_vars=groupby, value_vars=genes, var_name="gene", value_name="avg_expr")
    df_frac = frac_df.melt(id_vars=groupby, value_vars=genes, var_name="gene", value_name="frac_expr")

    # Merge
    df_plot = pd.merge(df_mean, df_frac, on=[groupby, "gene"])

    # Rename group column to 'group' for plotting
    df_plot = df_plot.rename(columns={groupby: "group"})

    return df_plot

