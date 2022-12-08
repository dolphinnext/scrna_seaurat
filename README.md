The standard Seurat scRNA downstream analysis pipeline  

For more information: 
https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

Parameters:
- minUMI: Minimum number of UMI counts for each cell (default: 1000)
- maxUMI: Maximum number of UMI counts for each cell (default: 10000)
- mitoRatio: Maximum mitochondrial contamination percentage for each cell (default: 5)
- varFeatures: The number of most variable genes (default: 2000, see Seurat tutorial)

Steps: 
- Filter cells with low and high UMI counts
- Normalize, and find most variable genes
- Scale UMI counts, find top PCs and cluster cells with resolution=0.4, 0.6, 0.8, 1.0, 1.4 (see Seurat tutorial)
- Project clusters to UMAP

Output: 
- Seurat object with resulting clusters
- h5ad file for Cellxgene and annotation
- UMAP of clusters
- UMAP of UMI counts
- UMI density plots of clusters
- Heatmap of marker genes (top 12)
- Table of markers genes (top 50)