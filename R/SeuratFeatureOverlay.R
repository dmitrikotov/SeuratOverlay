#' FeaturePlot with Cluster Overlays
#'
#' Generates a Seurat FeaturePlot with closed-loop boundaries around
#' clusters and custom repelled/static centroid labels. It filters out low-count
#' clusters and prunes border cells using a fast KNN homogeneity filter to prevent
#' overlapping boundary lines.
#'
#' @param seurat_obj A Seurat object.
#' @param gene_to_plot Character. Gene/feature to visualize.
#' @param reduction_name Character. Dimensional reduction to use (e.g., "umap", "tsne"). Default is "umap".
#' @param group_column Character. Metadata column containing cluster labels. Default is "seurat_clusters".
#' @param idents_to_plot Character vector. Specific identities within the group_column to outline and label. If NULL, all are plotted. Default is NULL.
#' @param group.by Character. Metadata column to group cells by in FeaturePlot. Default is NULL.
#' @param split.by Character. Metadata column to split the plot by in FeaturePlot. Default is NULL.
#' @param repel_labels Logical. Whether to repel labels away from centroids. Default is FALSE.
#' @param show_label_lines Logical. Whether to draw connector lines from repelled labels. Default is FALSE.
#' @param label_size Numeric. Font size for the cluster labels. Default is 4.
#' @param line_width Numeric. Stroke thickness of the boundaries. Default is 0.6.
#' @param line_type Character. Line type for boundaries (e.g., "dashed", "dotted", "solid"). Default is "dashed".
#' @param min_cells Numeric. Minimum cell count required to outline/label a cluster. Default is 50.
#' @param k_neighbors Numeric. Number of nearest neighbors to evaluate for border pruning. Default is 10.
#' @param neighbor_homogeneity Numeric. Value between 0-1. Required local purity to keep a cell for boundary drawing. Default is 0.85.
#'
#' @importFrom Seurat Embeddings FeaturePlot
#' @importFrom RANN nn2
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 geom_label geom_polygon theme_minimal labs aes unit
#' @importFrom dplyr %>% count filter group_by mutate slice summarize ungroup
#'
#' @return A ggplot object containing the styled FeaturePlot with overlays.
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' FeaturePlotWithOverlays(
#'   seurat_obj = pbmc_small,
#'   gene_to_plot = "CD8A",
#'   reduction_name = "tsne",
#'   group_column = "groups",
#'   min_cells = 30
#' )
#' }
FeaturePlotWithOverlays <- function(seurat_obj,
                                    gene_to_plot,
                                    reduction_name = "umap",
                                    group_column = "seurat_clusters",
                                    idents_to_plot = NULL,
                                    group.by = NULL,
                                    split.by = NULL,
                                    repel_labels = FALSE,
                                    show_label_lines = FALSE,
                                    label_size = 4,
                                    line_width = 0.6,
                                    line_type = "dashed",
                                    min_cells = 50,
                                    k_neighbors = 10,
                                    neighbor_homogeneity = 0.85) {

  # --- Step 1: Extract Coordinates and Metadata ---
  embed_coords <- as.data.frame(Seurat::Embeddings(seurat_obj, reduction = reduction_name))
  if (ncol(embed_coords) < 2) {
    stop("Selected dimensional reduction must have at least 2 dimensions.")
  }
  colnames(embed_coords)[1:2] <- c("Dim_1", "Dim_2")

  # Add the grouping variable from metadata
  if (!group_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Metadata column", group_column, "not found in Seurat object."))
  }
  embed_coords$Cluster <- seurat_obj[[group_column]][, 1]

  # --- Step 2: Validate split.by and group.by if provided ---
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(seurat_obj@meta.data)) {
      stop(paste("split.by column", split.by, "not found in Seurat object metadata."))
    }
  }
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(seurat_obj@meta.data)) {
      stop(paste("group.by column", group.by, "not found in Seurat object metadata."))
    }
  }

  # --- Step 3: Filter for Specific Idents if requested ---
  if (!is.null(idents_to_plot)) {
    invalid_idents <- setdiff(idents_to_plot, unique(embed_coords$Cluster))
    if (length(invalid_idents) > 0) {
      warning(paste("The following idents were not found in the metadata and will be ignored:",
                    paste(invalid_idents, collapse = ", ")))
    }
    embed_coords <- embed_coords %>% dplyr::filter(Cluster %in% idents_to_plot)
  }

  # --- Step 4: Filter out Small Clusters ---
  cluster_counts <- embed_coords %>%
    dplyr::count(Cluster) %>%
    dplyr::filter(n >= min_cells)

  if (nrow(cluster_counts) == 0) {
    stop(paste("No clusters have at least", min_cells, "cells. Try lowering min_cells."))
  }

  embed_coords <- embed_coords %>%
    dplyr::filter(Cluster %in% cluster_counts$Cluster)

  # --- Step 5: Fast KNN Neighborhood Homogeneity Filter ---
  if (nrow(embed_coords) > k_neighbors) {
    nn_results <- RANN::nn2(
      data = embed_coords[, c("Dim_1", "Dim_2")],
      k = k_neighbors + 1
    )

    # Extract neighbor indices, ignoring self (first column)
    neighbor_indices <- nn_results$nn.idx[, 2:(k_neighbors + 1)]

    # Map neighbor indices to corresponding cluster labels
    neighbor_clusters <- matrix(
      embed_coords$Cluster[neighbor_indices],
      nrow = nrow(embed_coords),
      ncol = k_neighbors
    )

    # Calculate homogeneity scores and filter out border/interface cells
    homogeneity_scores <- rowSums(neighbor_clusters == embed_coords$Cluster) / k_neighbors
    embed_coords$Homogeneity <- homogeneity_scores
    embed_coords_clean <- embed_coords %>% dplyr::filter(Homogeneity >= neighbor_homogeneity)
  } else {
    embed_coords_clean <- embed_coords
  }

  # --- Step 6: Compute Centroids for Labels ---
  centroids <- embed_coords %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(
      Dim_1 = mean(Dim_1),
      Dim_2 = mean(Dim_2)
    ) %>%
    dplyr::ungroup()

  # --- Step 7: Compute Convex Hulls for Boundary Outlines ---
  hull_data <- embed_coords_clean %>%
    dplyr::group_by(Cluster) %>%
    dplyr::mutate(
      Centroid_1 = mean(Dim_1),
      Centroid_2 = mean(Dim_2),
      Distance = sqrt((Dim_1 - Centroid_1)^2 + (Dim_2 - Centroid_2)^2)
    ) %>%
    dplyr::filter(Distance <= stats::quantile(Distance, 0.95)) %>% # Core 95%
    dplyr::slice(grDevices::chull(Dim_1, Dim_2)) %>%
    dplyr::ungroup()

  # --- Step 8: Generate Base FeaturePlot ---
  base_plot <- Seurat::FeaturePlot(
    object = seurat_obj,
    features = gene_to_plot,
    reduction = reduction_name,
    split.by = split.by,
    group.by = group.by
  )

  # --- Step 9: Configure Label Layer ---
  if (repel_labels) {
    label_layer <- ggrepel::geom_label_repel(
      data = centroids,
      ggplot2::aes(x = Dim_1, y = Dim_2, label = Cluster),
      fill = "white",
      color = "black",
      fontface = "bold",
      size = label_size,
      alpha = 0.85,
      label.padding = ggplot2::unit(0.2, "lines"),
      box.padding = 1.2,
      point.padding = 0, # Forces connector to attach directly to centroid
      force = 10,
      segment.color = if (show_label_lines) "grey30" else NA,
      segment.size = 0.5,
      min.segment.length = 0,
      arrow = NULL, # Removed arrow segment as requested
      inherit.aes = FALSE
    )
  } else {
    label_layer <- ggplot2::geom_label(
      data = centroids,
      ggplot2::aes(x = Dim_1, y = Dim_2, label = Cluster),
      fill = "white",
      color = "black",
      fontface = "bold",
      size = label_size,
      alpha = 0.85,
      label.padding = ggplot2::unit(0.2, "lines"),
      inherit.aes = FALSE
    )
  }

  # --- Step 10: Assemble and Return Final Plot ---
  final_plot <- base_plot +
    ggplot2::geom_polygon(
      data = hull_data,
      ggplot2::aes(x = Dim_1, y = Dim_2, group = Cluster),
      color = "grey30",
      fill = NA,
      linetype = line_type,
      linewidth = line_width,
      alpha = 0.75,
      inherit.aes = FALSE
    ) +
    label_layer +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Expression of", gene_to_plot),
      subtitle = paste0("Outlines show complete core areas (K=", k_neighbors, ", min homogeneity=", neighbor_homogeneity*100, "%)")
    )

  return(final_plot)
}
