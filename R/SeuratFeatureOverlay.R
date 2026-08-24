#' FeaturePlot with Cluster Overlays
#'
#' Generates a Seurat FeaturePlot with closed-loop boundaries around
#' clusters and custom repelled/static centroid labels. It filters out low-count
#' clusters and prunes border cells using a fast KNN homogeneity filter to prevent
#' overlapping boundary lines. Instead of a convex hull, it uses 2D kernel density
#' estimation to construct concave boundary lines that tightly follow the actual cell shapes.
#'
#' @param object A Seurat object.
#' @param features Character vector. Gene(s)/feature(s) to visualize.
#' @param reduction Character. Dimensional reduction to use (e.g., "umap", "tsne"). Default is "umap".
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
#' @param contour_threshold Numeric. Value between 0-1. The density threshold at which to draw the boundary lines, relative to the peak density of each cluster. Default is 0.05 (5 percent of peak density). Lower values wrap wider/looser; higher values hug tightly.
#' @param grid_size Numeric. Grid size (n x n) for the 2D kernel density estimation. Default is 100. Higher values yield smoother contours.
#' @param merge_threshold Numeric. Proximity threshold scaled as a fraction of the coordinate span (Dim 1 + Dim 2) below which separate segments of the same cluster are merged. Default is 0.08 (8 percent). Increase this value to merge segments that are further apart.
#' @param raster Logical. Whether to rasterize points in the Seurat FeaturePlot. Default is FALSE.
#' @param slot Character. Expression slot to pull data from (e.g., "data", "counts"). Default is "data".
#' @param assay Character. Specific assay to pull feature data from. Default is NULL.
#' @param largest_boundary_only Logical. If TRUE, retains and plots only the largest boundary polygon (by enclosed area) for each cluster with identical cluster names. Default is FALSE.
#'
#' @importFrom Seurat Embeddings FeaturePlot
#' @importFrom RANN nn2
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 geom_label geom_polygon theme_minimal labs aes unit theme guides guide_colorbar element_blank element_line element_text
#' @importFrom dplyr %>% count filter group_by mutate slice summarize ungroup group_modify slice_max
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @importFrom MASS kde2d
#' @importFrom grDevices contourLines chull
#' @importFrom mgcv in.out
#'
#' @return A ggplot object containing the styled FeaturePlot with overlays.
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' FeaturePlotWithOverlays(
#'   object = pbmc_small,
#'   features = c("CD8A", "GZMB"),
#'   reduction = "tsne",
#'   group_column = "groups",
#'   min_cells = 30
#' )
#' }
FeaturePlotWithOverlays <- function(object,
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
                                    largest_boundary_only = FALSE) {

  # Check if split.by is valid if provided (independent of group_column)
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(object@meta.data)) {
      stop(paste("split.by column", split.by, "not found in Seurat object metadata."))
    }
  }

  is_split <- !is.null(split.by)
  is_multi_gene <- length(features) > 1

  # Conditionally prepare overlay layers
  if (is.null(group_column)) {
    polygon_layer <- NULL
    label_layer   <- NULL
  } else {
    # --- Step 1: Extract Coordinates and Metadata ---
    embed_coords <- as.data.frame(Seurat::Embeddings(object, reduction = reduction))
    if (ncol(embed_coords) < 2) {
      stop("Selected dimensional reduction must have at least 2 dimensions.")
    }
    colnames(embed_coords)[1:2] <- c("Dim_1", "Dim_2")

    # Add the grouping variable from metadata
    if (!group_column %in% colnames(object@meta.data)) {
      stop(paste("Metadata column", group_column, "not found in Seurat object."))
    }
    embed_coords$Cluster <- object[[group_column]][, 1]

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

    # --- Step 5: Compute Concave Density-Based Outlines (Contours) ---
    # We compute the boundaries sequentially. For each cluster, we mask its density grid
    # to be 0 where other clusters are closer, AND where the grid points fall inside
    # any previously computed cluster boundaries. This guarantees that no cluster boundary
    # line can cross into or be inside another cluster's boundary.
    cluster_sizes <- table(embed_coords_clean$Cluster)
    unique_clusters <- names(sort(cluster_sizes, decreasing = FALSE)) # Smaller clusters first
    hull_list <- list()

    for (cl in unique_clusters) {
      cl_data <- embed_coords_clean[embed_coords_clean$Cluster == cl, ]
      if (nrow(cl_data) < 5) next

      # Calculate coordinate limits with 15% padding
      range_1 <- range(cl_data$Dim_1)
      range_2 <- range(cl_data$Dim_2)
      span_1 <- range_1[2] - range_1[1]
      span_2 <- range_2[2] - range_2[1]
      lims_1 <- c(range_1[1] - 0.15 * span_1, range_1[2] + 0.15 * span_1)
      lims_2 <- c(range_2[1] - 0.15 * span_2, range_2[2] + 0.15 * span_2)

      # Estimate 2D spatial density with padded limits
      dens <- MASS::kde2d(cl_data$Dim_1, cl_data$Dim_2, n = grid_size, lims = c(lims_1, lims_2))

      # Generate grid intersection coordinates
      grid_coords <- expand.grid(Dim_1 = dens$x, Dim_2 = dens$y)

      # Find the single closest real cell to each grid intersection across the clean dataset
      nn_grid <- RANN::nn2(
        data = embed_coords_clean[, c("Dim_1", "Dim_2")],
        query = grid_coords,
        k = 1
      )

      # Identify the closest cluster label for each grid coordinate
      nearest_cluster <- embed_coords_clean$Cluster[nn_grid$nn.idx[, 1]]

      # Mask density values to 0 along Voronoi decision boundaries with other clusters
      is_foreign <- nearest_cluster != cl
      dens$z[is_foreign] <- 0

      # ALSO mask density values to 0 if grid points fall inside any already computed boundaries.
      # This strictly prevents incoming cluster contours from nesting inside existing ones.
      if (length(hull_list) > 0) {
        for (prev_hull in hull_list) {
          for (p_id in unique(prev_hull$Piece)) {
            piece_data <- prev_hull[prev_hull$Piece == p_id, c("Dim_1", "Dim_2")]
            # Use mgcv::in.out for highly optimized C-compiled point-in-polygon test
            is_inside <- mgcv::in.out(as.matrix(piece_data), as.matrix(grid_coords))
            dens$z[is_inside] <- 0
          }
        }
      }

      # Find threshold value relative to peak density of the masked grid
      peak_density <- max(dens$z)
      level_threshold <- peak_density * contour_threshold

      if (level_threshold <= 0) next

      # Extract contour coordinates at threshold
      contours <- grDevices::contourLines(dens$x, dens$y, dens$z, levels = level_threshold)

      if (length(contours) == 0) next

      # Format individual contour polygons, ensuring closed paths.
      pieces_list <- lapply(seq_along(contours), function(i) {
        cx <- contours[[i]]$x
        cy <- contours[[i]]$y

        if (cx[1] != cx[length(cx)] || cy[1] != cy[length(cy)]) {
          cx <- c(cx, cx[1])
          cy <- c(cy, cy[1])
        }

        data.frame(
          Dim_1 = cx,
          Dim_2 = cy,
          Piece_ID = i
        )
      })

      num_pieces <- length(pieces_list)
      if (num_pieces == 0) next

      # Define spatial touching/merging threshold (fraction of total coordinate dimensions)
      touch_threshold <- merge_threshold * (span_1 + span_2)

      # Initialize each piece in its own group
      group_assignments <- seq_len(num_pieces)

      # Check pairwise proximity of segments to group touching/very close pieces together
      if (num_pieces > 1) {
        for (i in 1:(num_pieces - 1)) {
          for (j in (i + 1):num_pieces) {
            p_i <- pieces_list[[i]][, c("Dim_1", "Dim_2")]
            p_j <- pieces_list[[j]][, c("Dim_1", "Dim_2")]

            # Fast nearest neighbor search between contour boundaries
            nn_dist <- RANN::nn2(data = p_i, query = p_j, k = 1)
            min_d <- min(nn_dist$nn.dists)

            if (min_d < touch_threshold) {
              g_to_change <- group_assignments[j]
              group_assignments[group_assignments == g_to_change] <- group_assignments[i]
            }
          }
        }
      }

      # Assemble final polygons
      final_pieces <- list()
      unique_groups <- unique(group_assignments)

      for (g_idx in seq_along(unique_groups)) {
        g <- unique_groups[g_idx]
        member_indices <- which(group_assignments == g)

        # Combine all coordinates from pieces belonging to this group
        combined_pts <- do.call(rbind, lapply(member_indices, function(idx) pieces_list[[idx]]))

        if (length(member_indices) > 1) {
          # Joining boundaries by taking the convex hull of combined points of the touching pieces
          hull_idx <- grDevices::chull(combined_pts$Dim_1, combined_pts$Dim_2)
          joined_cx <- combined_pts$Dim_1[hull_idx]
          joined_cy <- combined_pts$Dim_2[hull_idx]

          # Ensure it forms a completely closed path
          joined_cx <- c(joined_cx, joined_cx[1])
          joined_cy <- c(joined_cy, joined_cy[1])

          final_pieces[[g_idx]] <- data.frame(
            Cluster = cl,
            Dim_1 = joined_cx,
            Dim_2 = joined_cy,
            Piece = paste0(cl, "_G", g_idx)
          )
        } else {
          # Keep original concave contour shape for non-touching islands
          final_pieces[[g_idx]] <- data.frame(
            Cluster = cl,
            Dim_1 = combined_pts$Dim_1,
            Dim_2 = combined_pts$Dim_2,
            Piece = paste0(cl, "_G", g_idx)
          )
        }
      }

      cluster_hulls <- do.call(rbind, final_pieces)
      if (!is.null(cluster_hulls) && nrow(cluster_hulls) > 0) {
        hull_list[[cl]] <- cluster_hulls
      }
    }

    # Combine all cluster hulls into a single data frame
    if (length(hull_list) > 0) {
      hull_data <- do.call(rbind, hull_list)
      rownames(hull_data) <- NULL
    } else {
      hull_data <- data.frame()
    }

    # Filter to only the single largest boundary piece per cluster if requested
    if (largest_boundary_only && nrow(hull_data) > 0) {
      poly_area <- function(x, y) {
        n <- length(x)
        if (n < 3) return(0)
        0.5 * abs(sum(x[1:(n - 1)] * y[2:n] - x[2:n] * y[1:(n - 1)]))
      }

      piece_areas <- hull_data %>%
        dplyr::group_by(Cluster, Piece) %>%
        dplyr::summarize(
          area = poly_area(Dim_1, Dim_2),
          .groups = "drop"
        ) %>%
        dplyr::group_by(Cluster) %>%
        dplyr::slice_max(area, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

      hull_data <- hull_data %>%
        dplyr::filter(Piece %in% piece_areas$Piece)
    }

    # --- Step 6: Compute Robust Centroids and Target Anchors for Labels ---
    if (nrow(hull_data) == 0) {
      polygon_layer <- NULL
      label_layer <- NULL
    } else {
      global_center <- c(mean(embed_coords$Dim_1), mean(embed_coords$Dim_2))

      centroids <- hull_data %>%
        dplyr::group_by(Cluster, Piece) %>%
        dplyr::group_modify(~ {
          poly_x <- .x$Dim_1
          poly_y <- .x$Dim_2
          center_x <- mean(poly_x)
          center_y <- mean(poly_y)

          # Get real cells belonging to this specific cluster
          cluster_cells <- embed_coords_clean %>% dplyr::filter(Cluster == .y$Cluster)

          if (nrow(cluster_cells) == 0) return(data.frame())

          # Find the robust medoid cell closest to the center of this specific boundary piece
          dists_to_center <- sqrt((cluster_cells$Dim_1 - center_x)^2 + (cluster_cells$Dim_2 - center_y)^2)
          best_cell_idx <- which.min(dists_to_center)[1]

          c_x <- cluster_cells$Dim_1[best_cell_idx]
          c_y <- cluster_cells$Dim_2[best_cell_idx]

          if (repel_labels) {
            # Find all other clusters' boundary points to repel away from
            other_hulls <- hull_data %>% dplyr::filter(Cluster != .y$Cluster)

            if (nrow(other_hulls) > 0) {
              # Fast pairwise calculation using nn2:
              # For each boundary coordinate of our current piece, calculate minimum distance to any other cluster border
              current_piece_coords <- data.frame(Dim_1 = poly_x, Dim_2 = poly_y)
              nn_dist_results <- RANN::nn2(
                data = other_hulls[, c("Dim_1", "Dim_2")],
                query = current_piece_coords,
                k = 1
              )

              # Find the index of our boundary point that maximizes the minimum distance to other clusters (argmax)
              best_idx <- which.max(nn_dist_results$nn.dists)[1]

              data.frame(
                Dim_1 = poly_x[best_idx],
                Dim_2 = poly_y[best_idx]
              )
            } else {
              # Fallback: Outward projection from the global coordinate center if no other clusters exist
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

              max_rad <- max(sqrt((cluster_cells$Dim_1 - c_x)^2 + (cluster_cells$Dim_2 - c_y)^2), na.rm = TRUE)
              est_label_x <- c_x + (max_rad * 1.5) * u_x
              est_label_y <- c_y + (max_rad * 1.5) * u_y

              dists <- sqrt((poly_x - est_label_x)^2 + (poly_y - est_label_y)^2)
              best_idx <- which.min(dists)[1]

              data.frame(
                Dim_1 = poly_x[best_idx],
                Dim_2 = poly_y[best_idx]
              )
            }
          } else {
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
          family = "Helvetica",
          size = label_size,
          alpha = 0.85,
          label.padding = ggplot2::unit(0.2, "lines"),
          box.padding = 1.2,
          point.padding = 0,
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
          family = "Helvetica",
          size = label_size,
          alpha = 0.85,
          label.padding = ggplot2::unit(0.2, "lines"),
          inherit.aes = FALSE
        )
      }

      # --- Step 8: Build Polygon Layer ---
      polygon_layer <- ggplot2::geom_polygon(
        data = hull_data,
        ggplot2::aes(x = Dim_1, y = Dim_2, group = Piece),
        color = "grey30",
        fill = NA,
        linetype = line_type,
        linewidth = line_width,
        alpha = 0.75,
        inherit.aes = FALSE
      )
    }
  }

  # --- Step 9: Generate Base FeaturePlot ---
  base_plot <- Seurat::FeaturePlot(
    object = object,
    features = features,
    reduction = reduction,
    split.by = split.by,
    order = order,
    raster = raster,
    slot = slot,
    assay = assay
  )

  # Standard guide settings to suppress legend titles for non-split plots
  null_guides <- ggplot2::guides(
    color = ggplot2::guide_colorbar(title = NULL)
  )

  # Theme modification for stripping grid lines, setting line widths, and using Helvetica with ticks.
  theme_clean_axes <- ggplot2::theme(
    text              = ggplot2::element_text(family = "Helvetica"),
    axis.text         = ggplot2::element_text(size = 12, color = "black", family = "Helvetica"),
    axis.text.y.right = ggplot2::element_blank(),
    axis.ticks.y.right= ggplot2::element_blank(),
    axis.title.y.right= ggplot2::element_blank(),
    axis.title        = ggplot2::element_text(family = "Helvetica", size = 14),
    plot.title        = ggplot2::element_text(family = "Helvetica"),
    legend.text       = ggplot2::element_text(family = "Helvetica", size = 11),
    legend.title      = ggplot2::element_text(family = "Helvetica", face = "bold.italic", size = 12),
    panel.grid.major  = ggplot2::element_blank(),
    panel.grid.minor  = ggplot2::element_blank(),
    axis.line         = ggplot2::element_line(color = "black", linewidth = 0.5),
    axis.line.y.right = ggplot2::element_blank(),
    axis.ticks        = ggplot2::element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = ggplot2::unit(0.15, "cm")
  )

  # --- Step 10: Assemble and Return Final Plot ---
  if (inherits(base_plot, "patchwork")) {
    if (is_split) {
      num_genes <- length(features)
      num_subplots <- length(base_plot)
      num_splits <- num_subplots / num_genes

      gene_row_plots <- list()

      for (g in seq_len(num_genes)) {
        indices <- ((g - 1) * num_splits + 1):(g * num_splits)

        # Define clean row guides displaying the specific feature/gene title in bold and italic
        clean_guides_g <- ggplot2::guides(
          color = ggplot2::guide_colorbar(title = features[g], title.position = "top")
        )

        # Apply overlays to the subplots of this specific gene
        for (i in indices) {
          base_plot[[i]] <- base_plot[[i]] +
            polygon_layer +
            label_layer +
            ggplot2::theme_minimal() +
            theme_clean_axes +
            ggplot2::theme(
              legend.position = "right",
              plot.title = ggplot2::element_text(face = "bold", size = 14, family = "Helvetica")
            ) +
            clean_guides_g
        }

        # Extract processed subplots for this gene
        gene_plots_list <- lapply(indices, function(idx) base_plot[[idx]])

        # Group panels, collect guides, and safely append annotation titles
        gene_row <- patchwork::wrap_plots(gene_plots_list, ncol = num_splits, guides = "collect") &
          ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(face = "bold", size = 14, family = "Helvetica")
          ) &
          theme_clean_axes &
          clean_guides_g

        gene_row_plots[[g]] <- gene_row +
          patchwork::plot_annotation(
            title = features[g],
            theme = ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold.italic", size = 16, family = "Helvetica")
            )
          )
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
          theme_clean_axes +
          ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(face = "bold.italic", size = 14, family = "Helvetica")
          ) +
          null_guides
      }

      final_plot <- base_plot +
        patchwork::plot_layout(guides = "keep") &
        ggplot2::theme(
          legend.position = "right",
          plot.title = ggplot2::element_text(face = "bold.italic", size = 14, family = "Helvetica")
        ) &
        theme_clean_axes &
        null_guides
    }
  } else {
    # Standard single-panel ggplot
    final_plot <- base_plot +
      polygon_layer +
      label_layer +
      ggplot2::theme_minimal() +
      theme_clean_axes +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold.italic", size = 16, family = "Helvetica")
      ) +
      null_guides +
      ggplot2::labs(title = features)
  }

  return(final_plot)
}
