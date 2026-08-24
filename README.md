# SeuratOverlay
R package to overlay cluster metadata on top of the Seurat FeaturePlot function

# Installing SeuratOverlay
To install SeuratOverlay and the required dependencies from Github run the following code:

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
  object,
  features,
  reduction = "umap",
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
  neighbor_homogeneity = 0.85,
  contour_threshold = 0.05,
  grid_size = 100,
  merge_threshold = 0.08,
  raster = FALSE,
  slot = "data",
  assay = NULL,
  largest_boundary_only = FALSE
)

object: A Seurat object.

features: Gene/feature to visualize.

reduction: Dimensional reduction to use (e.g., "umap", "tsne"). Default is "umap".

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

contour_threshold: Value between 0-1. The density threshold at which to draw the boundary lines, relative to the peak density of each cluster. Default is 0.05 (5 percent of peak density). Lower values wrap wider/looser; higher values hug tightly.

grid_size: Grid size (n x n) for the 2D kernel density estimation. Default is 100. Higher values yield smoother contours.

merge_threshold: Proximity threshold scaled as a fraction of the coordinate span (Dim 1 + Dim 2) below which separate segments of the same cluster are merged. Default is 0.08 (8 percent). Increase this value to merge segments that are further apart.

raster: Whether to rasterize points in the Seurat FeaturePlot. Default is FALSE.

slot: Expression slot to pull data from (e.g., "data", "counts"). Default is "data".

assay: Specific assay to pull feature data from. Default is NULL.

largest_boundary_only: If TRUE, retains and plots only the largest boundary polygon (by enclosed area) for each cluster with identical cluster names. Default is FALSE.
```
