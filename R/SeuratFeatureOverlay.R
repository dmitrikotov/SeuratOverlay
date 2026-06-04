#' FeaturePlot with Cluster Overlays
#'
#' Generates a Seurat FeaturePlot with closed-loop boundaries around
#' clusters and custom repelled/static centroid labels. It filters out low-count
#' clusters and prunes border cells using a fast KNN homogeneity filter to prevent
#' overlapping boundary lines.
#'
#' @param seurat_obj A Seurat object.
#' @param gene_to_plot Character vector. Gene(s)/feature(s) to visualize.
#' @param reduction_name Character. Dimensional reduction to use (e.g., "umap", "tsne"). Default is "umap".
#' @param group_column Character. Metadata column containing cluster labels. If NULL, overlays are skipped. Default is "seurat_clusters".
#' @param idents_to_plot Character vector. Specific identities within the group_column to outline and label. If NULL, all are plotted. Default is NULL.
#' @param split.by Character. Metadata column to split the plot by in FeaturePlot. Default is NULL.
#' @param order Logical. If TRUE, plot cells in order of expression level so that cells expressing the feature are plotted on top. Default is TRUE.
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
#' @importFrom ggplot2 geom_label geom_polygon theme_minimal labs aes unit theme guides guide_colorbar
#' @importFrom dplyr %>% count filter group_by mutate slice summarize ungroup group_modify
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
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
#'   gene_to_plot = c("CD8A", "GZMB"),
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
                                    split.by = NULL,
                                    order = TRUE,
                                    repel_labels = FALSE,
                                    show_label_lines = FALSE,
                                    label_size = 4,
                                    line_width = 0.6,
                                    line_type = "dashed",
                                    min_cells = 50,
                                    k_neighbors = 10,
                                    neighbor_homogeneity = 0.85) {

  # Check if split.by is valid if provided (independent of group_column)
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(seurat_obj@meta.data)) {
      stop(paste("split.by column", split.by, "not found in Seurat object metadata."))
    }
  }

  is_split <- !is.null(split.by)
  is_multi_gene <- length(gene_to_plot) > 1

  # Conditionally prepare overlay layers
  if (is.null(group_column)) {
    polygon_layer <- NULL
    label_layer   <- NULL
  } else {
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

    # --- Step 2: Filter for Specific Idents if requested ---
    if (!is.null(idents_to_plot)) {
      invalid_idents <- setdiff(idents_to_plot, unique(embed_coords$Cluster))
      if (length(invalid_idents) > 0) {
        warning(paste("The following idents were not found in the metadata and will be ignored:",
                      paste(invalid_idents, collapse = ", ")))
      }
      embed_coords <- embed_coords %>% dplyr::filter(Cluster %in% idents_to_plot)
    }

    # --- Step 3: Filter out Small Clusters ---
    cluster_counts <- embed_coords %>%
      dplyr::count(Cluster) %>%
      dplyr::filter(n >= min_cells)

    if (nrow(cluster_counts) == 0) {
      stop(paste("No clusters have at least", min_cells, "cells. Try lowering min_cells."))
    }

    embed_coords <- embed_coords %>%
      dplyr::filter(Cluster %in% cluster_counts$Cluster)

    # --- Step 4: Fast KNN Neighborhood Homogeneity Filter ---
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

    # --- Step 5: Compute Convex Hulls for Boundary Outlines ---
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

    # --- Step 6: Compute Centroids and Target Anchors for Labels ---
    # We calculate the global center of coordinates to determine outward projection direction.
    global_center <- c(mean(embed_coords$Dim_1), mean(embed_coords$Dim_2))

    centroids <- hull_data %>%
      dplyr::group_by(Cluster) %>%
      dplyr::group_modify(~ {
        # Raw cluster centroid
        c_x <- mean(.x$Dim_1)
        c_y <- mean(.x$Dim_2)

        if (repel_labels) {
          # Calculate outward vector from global coordinate center
          v_x <- c_x - global_center[1]
          v_y <- c_y - global_center[2]
          v_len <- sqrt(v_x^2 + v_y^2)

          if (v_len > 0) {
            u_x <- v_x / v_len
            u_y <- v_y / v_len
          } else {
            u_x <- 0
            u_y <- 1
          }

          # Project a point outward from centroid along the unit vector
          max_rad <- max(sqrt((.x$Dim_1 - c_x)^2 + (.x$Dim_2 - c_y)^2))
          est_label_x <- c_x + (max_rad * 1.5) * u_x
          est_label_y <- c_y + (max_rad * 1.5) * u_y

          # Identify the point on the convex hull closest to this outward projected target
          dists <- sqrt((.x$Dim_1 - est_label_x)^2 + (.x$Dim_2 - est_label_y)^2)
          best_idx <- which.min(dists)

          # Use the closest hull boundary point as the anchor coordinate for repelled connector lines
          data.frame(
            Dim_1 = .x$Dim_1[best_idx],
            Dim_2 = .x$Dim_2[best_idx]
          )
        } else {
          # Keep coordinate centered exactly at raw centroid for non-repelled labels
          data.frame(
            Dim_1 = c_x,
            Dim_2 = c_y
          )
        }
      }) %>%
      dplyr::ungroup()

    # --- Step 7: Configure Label Layer ---
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
        point.padding = 0, # Forces connector to attach directly to boundary anchor
        force = 10,
        segment.color = if (show_label_lines) "grey30" else NA,
        segment.size = 0.5,
        min.segment.length = 0,
        arrow = NULL,
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

    # --- Step 8: Build Polygon Layer ---
    polygon_layer <- ggplot2::geom_polygon(
      data = hull_data,
      ggplot2::aes(x = Dim_1, y = Dim_2, group = Cluster),
      color = "grey30",
      fill = NA,
      linetype = line_type,
      linewidth = line_width,
      alpha = 0.75,
      inherit.aes = FALSE
    )
  }

  # --- Step 9: Generate Base FeaturePlot ---
  base_plot <- Seurat::FeaturePlot(
    object = seurat_obj,
    features = gene_to_plot,
    reduction = reduction_name,
    split.by = split.by,
    order = order
  )

  # Set up a target guide setup to kill titles using only the standardized 'color' aesthetic
  clean_guides <- ggplot2::guides(
    color = ggplot2::guide_colorbar(title = NULL)
  )

  # --- Step 10: Assemble and Return Final Plot ---
  # Apply overlays conditionally based on layout type (NULL layers are safely ignored by ggplot2)
  if (inherits(base_plot, "patchwork")) {
    if (is_split) {
      num_genes <- length(gene_to_plot)
      num_subplots <- length(base_plot)
      num_splits <- num_subplots / num_genes

      gene_row_plots <- list()

      for (g in seq_len(num_genes)) {
        indices <- ((g - 1) * num_splits + 1):(g * num_splits)

        # Apply overlays to the subplots of this specific gene
        for (i in indices) {
          base_plot[[i]] <- base_plot[[i]] +
            polygon_layer +
            label_layer +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "right") +
            clean_guides
        }

        # Extract processed subplots for this gene
        gene_plots_list <- lapply(indices, function(idx) base_plot[[idx]])

        # Group panels, collect guides, and safely append annotation titles
        gene_row <- patchwork::wrap_plots(gene_plots_list, guides = "collect") &
          ggplot2::theme(legend.position = "right") &
          clean_guides

        gene_row_plots[[g]] <- gene_row + patchwork::plot_annotation(title = gene_to_plot[g])
      }

      # Stack all clean rows vertically
      final_plot <- patchwork::wrap_plots(gene_row_plots, ncol = 1)

    } else {
      # Multi-gene plot with NO split.by (regular panel grid)
      for (i in seq_along(base_plot)) {
        base_plot[[i]] <- base_plot[[i]] +
          polygon_layer +
          label_layer +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "right") +
          clean_guides
      }

      # Keep legends corresponding to each individual gene subplot
      final_plot <- base_plot +
        patchwork::plot_layout(guides = "keep") &
        ggplot2::theme(legend.position = "right") &
        clean_guides
    }
  } else {
    # Standard single-panel ggplot (single gene, no split.by)
    final_plot <- base_plot +
      polygon_layer +
      label_layer +
      ggplot2::theme_minimal() +
      clean_guides +
      ggplot2::labs(title = gene_to_plot)
  }

  return(final_plot)
}
