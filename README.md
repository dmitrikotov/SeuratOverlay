# SeuratOverlay
R package to overlay cluster metadata on top of the Seurat FeaturePlot function

# Installing SeuratOverlay
To install SeuratOverlay and the required depndencies from github run the following code:

```R
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}

devtools::install_github("dmitrikotov/SeuratOverlay")
```

The main function within the SeuratOverlay package is FeaturePlotWithOverlays. This is a layer that runs on top of the Seurat FeaturePlot function allowing simultaneous plotting of gene expression data and cluster metadata.
Key inputs and features of the function are as follows:

```R
FeaturePlotWithOverlays(
  seurat_obj,
  gene_to_plot,
  reduction_name = "umap",
  group_column = "seurat_clusters",
  idents_to_plot = NULL,
  split.by = NULL,
  order = TRUE,
  repel_labels = FALSE,
  show_label_lines = FALSE,
  label_size = 4,
  line_width = 0.6,
  line_type = "dashed",
  min_cells = 50,
  k_neighbors = 10,
  neighbor_homogeneity = 0.85
)

seurat_obj: A Seurat object.

gene_to_plot: Gene/feature to visualize.

reduction_name: Dimensional reduction to use (e.g., "umap", "tsne"). Default is "umap".

group_column: Metadata column containing cluster labels. If NULL, overlays are skipped. Default is "seurat_clusters".

idents_to_plot: Specific identities within the group_column to outline and label. If NULL, all are plotted. Default is NULL.

split.by: Metadata column to split the plot by in FeaturePlot. Default is NULL.

order: If TRUE, plot cells in order of expression level so that cells expressing the feature are plotted on top. Default is TRUE.

repel_labels: Whether to repel labels away from centroids. Default is FALSE.

show_label_lines: Whether to draw connector lines from repelled labels. Default is FALSE.

label_size: Font size for the cluster labels. Default is 4.

line_width: Stroke thickness of the boundaries. Default is 0.6.

line_type: Line type for boundaries (e.g., "dashed", "dotted", "solid"). Default is "dashed".

min_cells: Minimum cell count required to outline/label a cluster. Default is 50.

k_neighbors: Number of nearest neighbors to evaluate for border pruning. Default is 10.

neighbor_homogeneity: Value between 0-1. Required local purity to keep a cell for boundary drawing. Default is 0.85.
```
