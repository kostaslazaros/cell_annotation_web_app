import scanpy as sc
import numpy as np
import pandas as pd
import hdf5plugin
import scipy as sci
import matplotlib.pyplot as plt
import seaborn as sns


def adata_preprocessor(
        adata_path, 
        min_genes=100, 
        min_cells=3, 
        n_genes_min=1000, 
        n_genes_max=10000, 
        n_counts_max=30000, 
        pc_mito=20, 
        pc_rib=25,
        debug=True
):
    """
    min_genes    : minimum genes that should be expressed in each cell
    min_cells    : minimum number of cells expected to express each gene
    n_genes_min  : lower cut-off of gene number for filtering cells
    n_genes_max  : upper cut-off of gene number for filtering cells
    n_counts_max : upper cut-off of counts for filtering cells
    pc_mito      : percentage of mitochondrial gene expression (cut-off)
    pc_rib       : percentage of ribosomal gene expression (cut-off)
    """
    
    def prn(txt):
        if debug:
            print(txt)
    
    fadata = sc.read_h5ad(adata_path)
    prn(f"Data shape before preprocessing: {fadata.shape}")
    
    # Pre-filtering
    sc.pp.filter_cells(fadata, min_genes=min_genes)  # Equivalent to min.features in Seurat.
    prn(f"Filtering cells with number of genes < 100: {fadata.shape}")
    
    sc.pp.filter_genes(fadata, min_cells=min_cells)  # Equivalent to min.cells in Seurat.
    prn(f"Filtering genes expressed in < 3 cells: {fadata.shape}")

    # Calculate the percentage of mitochondrial genes.
    mito_genes = fadata.var_names.str.startswith(tuple(['MT-', 'mt-', 'MT.', "mt."]))
    fadata.obs['prc_mt'] = (fadata[:, mito_genes].X.sum(axis=1) / fadata.X.sum(axis=1)) * 100
    prn("Mitochondrial gene percentage calculated and annotated in the prc_mt observation")

    # Calculate the percentage of ribosomal genes.
    ribo_genes = fadata.var_names.str.startswith('RPS')
    fadata.obs['prc_rb'] = (fadata[:, ribo_genes].X.sum(axis=1) / fadata.X.sum(axis=1)) * 100
    prn("Ribosomal gene percentage calculated and annotated in the prc_rb observation")

    # Calculate number of genes and counts for each cell.
    fadata.obs['n_genes'] = (fadata.X > 0).sum(axis=1)
    prn("Calculate number of genes with non-zero counts")
    
    fadata.obs['n_counts'] = fadata.X.sum(axis=1)
    prn("Calculate total number of counts for each cell")

    # Subsetting the data based on the calculated values.
    fadata = fadata[fadata.obs['n_genes'] > n_genes_min, :]
    prn(f"Filter cells with too few genes detected: {fadata.shape}")
    
    fadata = fadata[fadata.obs['n_genes'] < n_genes_max, :]
    prn(f"Filter cells with too many genes detected: {fadata.shape}")
    
    fadata = fadata[fadata.obs['n_counts'] < n_counts_max, :]
    prn(f"Filter cells with too many counts detected: {fadata.shape}")
    
    fadata = fadata[fadata.obs['prc_mt'] < pc_mito, :]
    prn(f"Filter cells with too many mitochondrial genes expressed: {fadata.shape}")
    
    fadata = fadata[fadata.obs['prc_rb'] < pc_rib, :]
    prn(f"Filter cells with too many ribosomal genes expressed: {fadata.shape}")
    
    return fadata


def sanitize_cellmarker_df(cellmarker_df, tissue_label, cols, drop_col_no=1):
    filtered_cellmarker_df = cellmarker_df[cellmarker_df[tissue_label[0]] == tissue_label[1]]
    filtered_cellmarker_df = filtered_cellmarker_df.loc[:, cols]
    filtered_cellmarker_df = filtered_cellmarker_df.dropna(subset=[cols[drop_col_no]])
    return filtered_cellmarker_df


def create_marker_dict(df_cellmarker, cell_name="Cell_type", marker="Symbol"):
    cell_markers = {}
    for i in list(set(df_cellmarker[cell_name].values.tolist())):
        genes = df_cellmarker[df_cellmarker[cell_name]==i][marker].values.tolist()
        cell_markers[i] = genes
    return cell_markers


def create_melted_df(score_df, ctype_lst):
    # Step 1: Creating separate DataFrames for each cluster
    cluster_dfs = [score_df[score_df['cluster'] == str(i)] for i in range(len(score_df["cluster"].value_counts()))]

    # Step 2: Adding a 'dataset' column
    dfs = [df.assign(dataset=f'dataset_{i}') for i, df in enumerate(cluster_dfs)]

    # Checking the assignment of the 'dataset' column

    # Step 3: Concatenating DataFrames into one
    final_df = pd.concat(dfs, ignore_index=True)

    # Confirming that the 'dataset' column is present in the concatenated DataFrame


    # Step 4: Melting the DataFrame for seaborn plotting
    final_df_melted = final_df.melt(id_vars=['dataset', 'cluster'], 
                                    value_vars=ctype_lst,
                                    var_name='cell_type',
                                    value_name='score')

    # Ensuring the melted DataFrame looks as expected
    return final_df_melted
    
    
def save_adata(adata_object, output_path):
    adata_object.write_h5ad(
        output_path,
        compression=hdf5plugin.FILTERS["zstd"]
    )
