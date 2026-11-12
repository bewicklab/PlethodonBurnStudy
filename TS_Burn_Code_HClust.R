
#===============================================================================
# MICROBIAL COMMUNITY HIERARCHICAL CLUSTERING
#
# Distance metrics:
#   1. Weighted UniFrac
#   2. Unweighted UniFrac
#   3. Jaccard
#   4. Bray-Curtis
#
# Clustering: Complete linkage
#
# Outputs:
#   1. Standalone weighted UniFrac dendrogram (unmodified)
#   2. Three vertically stacked dendrograms:
#      A. Jaccard
#      B. Bray-Curtis
#      C. Unweighted UniFrac
#
# In the three-panel figure:
#   - All sample tips remain aligned.
#   - The terminal region is compressed.
#   - Internal clustering topology is preserved.
#===============================================================================


#-------------------------------------------------------------------------------
# 1. PACKAGES AND SETTINGS
#-------------------------------------------------------------------------------

library(qiime2R)
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggtree)
library(ape)
library(patchwork)


output_dir <- "Composition_Analysis_Results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

set.seed(123)

# Terminal region compression settings
terminal_start_fraction <- 0.80
terminal_compression <- 0.15

# Site colors
site_colors <- c(
  "Control_1" = "#228A22",
  "Control_2" = "#00FE00",
  "Burn" = "red"
)

# Year shapes
year_shapes <- c(
  "2024" = 16,
  "2025" = 17
)


#-------------------------------------------------------------------------------
# 2. LOAD AND CLEAN METADATA
#-------------------------------------------------------------------------------

meta <- read.csv(
  "Coweeta_TEWL_Seasons_Sites.csv",
  stringsAsFactors = FALSE
)

meta <- meta[
  !is.na(meta$QIIME_ID_aligned) &
    meta$QIIME_ID_aligned != "",
]

meta <- meta[
  !duplicated(meta$QIIME_ID_aligned),
]

rownames(meta) <- tolower(
  gsub("-", ".", meta$QIIME_ID_aligned)
)


#-------------------------------------------------------------------------------
# 3. IMPORT QIIME2 DATA
#-------------------------------------------------------------------------------

Z_physeq <- qza_to_phyloseq(
  features = "table.qza",
  taxonomy = "taxonomy.qza",
  tree = "tree.qza"
)

sample_names(Z_physeq) <- tolower(
  gsub("-", ".", sample_names(Z_physeq))
)

common_samples <- intersect(
  sample_names(Z_physeq),
  rownames(meta)
)

Z_physeq <- prune_samples(
  common_samples,
  Z_physeq
)

meta <- meta[
  common_samples,
  ,
  drop = FALSE
]

sample_data(Z_physeq) <- sample_data(meta)

stopifnot(
  all(
    sample_names(Z_physeq) ==
      rownames(sample_data(Z_physeq))
  )
)


#-------------------------------------------------------------------------------
# 4. REMOVE LOW-ABUNDANCE TAXA AND SUBTRACT BLANK CONTROLS
#-------------------------------------------------------------------------------

Z_physeq <- prune_taxa(
  taxa_sums(Z_physeq) > 10,
  Z_physeq
)

control_labels <- "Blank"

control_samples <- sample_names(Z_physeq)[
  sample_data(Z_physeq)$Site %in% control_labels
]

if (length(control_samples) == 0) {
  stop("No blank control samples detected.")
}

otu_mat <- as(
  otu_table(Z_physeq),
  "matrix"
)

if (!taxa_are_rows(Z_physeq)) {
  otu_mat <- t(otu_mat)
}

control_mat <- otu_mat[
  ,
  control_samples,
  drop = FALSE
]

max_control_reads <- apply(
  control_mat,
  1,
  max,
  na.rm = TRUE
)

max_control_reads[
  !is.finite(max_control_reads)
] <- 0

corrected_mat <- sweep(
  otu_mat,
  1,
  max_control_reads,
  FUN = "-"
)

corrected_mat[corrected_mat < 0] <- 0

otu_table(Z_physeq) <- otu_table(
  corrected_mat,
  taxa_are_rows = TRUE
)

Z_physeq <- prune_samples(
  !(sample_data(Z_physeq)$Site %in% control_labels),
  Z_physeq
)

Z_physeq <- prune_taxa(
  taxa_sums(Z_physeq) > 0,
  Z_physeq
)


#-------------------------------------------------------------------------------
# 5. REMOVE UNWANTED TAXA
#-------------------------------------------------------------------------------

Z_physeq <- subset_taxa(
  Z_physeq,
  !(Kingdom %in% "d__Eukaryota") &
    !(Phylum %in% c(
      "Chytridiomycota",
      "Ascomycota",
      "Chlorophyta",
      "Apicomplexa",
      "Ciliophora",
      "Phragmoplastophyta",
      "Zoopagomycota"
    )) &
    Family != "Chloroplast" &
    Genus != "Mitochondria"
)

Z_physeq <- prune_taxa(
  taxa_sums(Z_physeq) > 0,
  Z_physeq
)


#-------------------------------------------------------------------------------
# 6. REMOVE UNWANTED SAMPLES
#-------------------------------------------------------------------------------

excluded_sites <- c(
  "Control",
  "2750",
  "2850",
  "EF",
  "PN",
  "Positive_control"
)

keep_samples <- !(
  sample_data(Z_physeq)$Site %in% excluded_sites
)

excluded_titles <- paste0(
  "TS.31.",
  5:11
)

keep_samples <- keep_samples &
  !(
    sample_data(Z_physeq)$title %in% excluded_titles
  )

excluded_seasons <- c(
  "Summer_2024",
  "Summer_2025",
  "Fall_2024"
)

keep_samples <- keep_samples &
  !(
    sample_data(Z_physeq)$Season %in% excluded_seasons
  )

keep_samples <- keep_samples &
  !is.na(
    sample_data(Z_physeq)$lat_7_Collection_locatio
  )

Z_physeq <- prune_samples(
  keep_samples,
  Z_physeq
)

Z_physeq <- prune_taxa(
  taxa_sums(Z_physeq) > 0,
  Z_physeq
)


#-------------------------------------------------------------------------------
# 7. RAREFACTION
#-------------------------------------------------------------------------------

rZOTU <- rarefy_even_depth(
  Z_physeq,
  sample.size = 12000,
  rngseed = 123,
  replace = FALSE
)


#-------------------------------------------------------------------------------
# 8. CORRECT SITE NAMES AND PREPARE METADATA
#-------------------------------------------------------------------------------

site <- as.character(
  sample_data(rZOTU)$Site
)

site[site == "31_burn"] <- "Burn"
site[site == "32_no_burn"] <- "Control_2"
site[site == "2950"] <- "Control_1"

sample_data(rZOTU)$Site <- factor(
  site,
  levels = c(
    "Control_1",
    "Control_2",
    "Burn"
  )
)

year <- sub(
  ".*_(20[0-9]{2})$",
  "\\1",
  as.character(sample_data(rZOTU)$Season)
)

sample_data(rZOTU)$Year <- year

tip_data <- data.frame(
  Site = as.character(
    sample_data(rZOTU)$Site
  ),
  Year = year,
  row.names = sample_names(rZOTU)
)

if (anyNA(tip_data$Site) || anyNA(tip_data$Year)) {
  stop("Missing site or year information.")
}

cat("\nSample sizes:\n")
print(table(tip_data$Site, tip_data$Year))


#-------------------------------------------------------------------------------
# 9. CALCULATE DISTANCE MATRICES
#-------------------------------------------------------------------------------

calculate_distance <- function(metric) {
  
  if (metric == "jaccard") {
    
    comm <- as(
      otu_table(rZOTU),
      "matrix"
    )
    
    if (taxa_are_rows(rZOTU)) {
      comm <- t(comm)
    }
    
    d <- vegan::vegdist(
      comm > 0,
      method = "jaccard",
      binary = TRUE
    )
    
  } else if (metric == "bray") {
    
    comm <- as(
      otu_table(rZOTU),
      "matrix"
    )
    
    if (taxa_are_rows(rZOTU)) {
      comm <- t(comm)
    }
    
    d <- vegan::vegdist(
      comm,
      method = "bray"
    )
    
  } else if (metric == "unifrac") {
    
    d <- phyloseq::distance(
      rZOTU,
      method = "unifrac",
      weighted = FALSE
    )
    
  } else if (metric == "wunifrac") {
    
    d <- phyloseq::distance(
      rZOTU,
      method = "unifrac",
      weighted = TRUE
    )
    
  } else {
    
    stop("Unknown distance metric: ", metric)
    
  }
  
  dmat <- as.matrix(d)
  
  sample_ids <- sample_names(rZOTU)
  
  if (
    !all(sample_ids %in% rownames(dmat)) ||
    !all(sample_ids %in% colnames(dmat))
  ) {
    stop("Distance matrix sample mismatch: ", metric)
  }
  
  dmat <- dmat[
    sample_ids,
    sample_ids,
    drop = FALSE
  ]
  
  if (any(!is.finite(dmat))) {
    stop("Non-finite distances detected: ", metric)
  }
  
  as.dist(dmat)
}


#-------------------------------------------------------------------------------
# 10. COMPRESS TERMINAL REGION WHILE KEEPING ALL TIPS ALIGNED
#-------------------------------------------------------------------------------

compress_terminal_region <- function(
    tree,
    start_fraction = 0.80,
    compression = 0.15) {
  
  if (
    start_fraction <= 0 ||
    start_fraction >= 1 ||
    compression <= 0 ||
    compression > 1
  ) {
    stop("Invalid terminal compression settings.")
  }
  
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths to compress.")
  }
  
  # Calculate original distance from root to every node
  original_depth <- ape::node.depth.edgelength(tree)
  
  # Maximum root-to-tip depth
  max_depth <- max(
    original_depth[seq_len(length(tree$tip.label))]
  )
  
  # Position where compression begins
  cutoff <- max_depth * start_fraction
  
  # Compress the portion beyond the cutoff
  new_depth <- ifelse(
    original_depth <= cutoff,
    original_depth,
    cutoff + (original_depth - cutoff) * compression
  )
  
  # Recalculate every edge length using transformed node depths
  parent_nodes <- tree$edge[, 1]
  child_nodes <- tree$edge[, 2]
  
  new_edge_lengths <- (
    new_depth[child_nodes] -
      new_depth[parent_nodes]
  )
  
  # Numerical protection against tiny negative values
  new_edge_lengths[
    new_edge_lengths < 0 &
      new_edge_lengths > -1e-10
  ] <- 0
  
  if (any(new_edge_lengths < 0)) {
    stop("Negative branch lengths generated during compression.")
  }
  
  tree$edge.length <- new_edge_lengths
  
  # Verify that all tips remain at the same endpoint
  compressed_depth <- ape::node.depth.edgelength(tree)
  
  compressed_tip_depth <- compressed_depth[
    seq_len(length(tree$tip.label))
  ]
  
  if (
    diff(range(compressed_tip_depth)) > 1e-8
  ) {
    stop("Compressed tree tips are not aligned.")
  }
  
  tree
}


#-------------------------------------------------------------------------------
# 11. FUNCTION TO GENERATE A DENDROGRAM
#-------------------------------------------------------------------------------

make_dendrogram <- function(
    metric,
    metric_label,
    compress_terminal = FALSE) {
  
  cat("\nCalculating:", metric_label, "\n")
  
  # Calculate distance matrix
  dist_mat <- calculate_distance(metric)
  
  # Complete-linkage hierarchical clustering
  hclust_result <- hclust(
    dist_mat,
    method = "complete"
  )
  
  # Convert to phylogenetic tree
  tree <- ape::as.phylo(hclust_result)
  
  # Compress terminal region only for selected panels
  if (compress_terminal) {
    
    tree <- compress_terminal_region(
      tree,
      start_fraction = terminal_start_fraction,
      compression = terminal_compression
    )
    
    cat(
      "Terminal region compressed:",
      terminal_compression * 100,
      "% of original length beyond",
      terminal_start_fraction * 100,
      "% of tree depth.\n"
    )
  }
  
  # Create ggtree object
  p <- ggtree::ggtree(tree)
  
  # Match tips to sample metadata
  tip_indices <- match(
    p$data$label,
    rownames(tip_data)
  )
  
  if (anyNA(tip_indices[p$data$isTip])) {
    stop(
      "Could not match all tree tips to metadata for ",
      metric
    )
  }
  
  p$data$Site <- tip_data$Site[tip_indices]
  p$data$Year <- tip_data$Year[tip_indices]
  
  # Plot
  p <- p +
    geom_tippoint(
      aes(
        color = Site,
        shape = Year
      ),
      size = 3
    ) +
    scale_color_manual(
      values = site_colors,
      breaks = names(site_colors),
      labels = c(
        "Control 1",
        "Control 2",
        "Burn"
      ),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = year_shapes,
      breaks = names(year_shapes),
      drop = FALSE
    ) +
    labs(
      title = metric_label,
      color = "Site",
      shape = "Year"
    ) +
    coord_flip() +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_tree2() +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),
      legend.position = "right",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  p
}


#-------------------------------------------------------------------------------
# 12. GENERATE ALL FOUR DENDROGRAMS
#-------------------------------------------------------------------------------

# Weighted UniFrac: original branch lengths
p_wunifrac <- make_dendrogram(
  metric = "wunifrac",
  metric_label = "Weighted UniFrac",
  compress_terminal = FALSE
)

# Three-panel figure: compressed terminal regions
p_jaccard <- make_dendrogram(
  metric = "jaccard",
  metric_label = "Jaccard",
  compress_terminal = TRUE
)

p_bray <- make_dendrogram(
  metric = "bray",
  metric_label = "Bray-Curtis",
  compress_terminal = TRUE
)

p_unifrac <- make_dendrogram(
  metric = "unifrac",
  metric_label = "Unweighted UniFrac",
  compress_terminal = TRUE
)


#-------------------------------------------------------------------------------
# 13. SAVE WEIGHTED UNIFRAC FIGURE
#-------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "microbial_dendrogram_weighted_unifrac.pdf"
  ),
  plot = p_wunifrac,
  width = 12,
  height = 8,
  units = "in"
)


#-------------------------------------------------------------------------------
# 14. CREATE THREE VERTICALLY STACKED PANELS
#-------------------------------------------------------------------------------

three_panel_figure <- (
  p_jaccard /
    p_bray /
    p_unifrac
) +
  plot_layout(
    ncol = 1,
    guides = "collect",
    heights = c(1, 1, 1)
  ) &
  theme(
    legend.position = "right"
  )

print(three_panel_figure)


#-------------------------------------------------------------------------------
# 15. SAVE THREE-PANEL PDF
#-------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "microbial_dendrogram_three_metrics_vertical.pdf"
  ),
  plot = three_panel_figure,
  width = 12,
  height = 20,
  units = "in",
  limitsize = FALSE
)


#-------------------------------------------------------------------------------
# 16. SAVE SAMPLE METADATA
#-------------------------------------------------------------------------------

write.csv(
  tip_data,
  file = file.path(
    output_dir,
    "Dendrogram_sample_metadata.csv"
  ),
  row.names = TRUE
)


#-------------------------------------------------------------------------------
# 17. FINAL SUMMARY
#-------------------------------------------------------------------------------

cat("\n========================================\n")
cat("HIERARCHICAL CLUSTERING COMPLETE\n")
cat("========================================\n")

cat("\nSamples included:", nsamples(rZOTU), "\n")

cat("\nOutputs saved in:", output_dir, "\n")

cat("\n1. microbial_dendrogram_weighted_unifrac.pdf\n")
cat("2. microbial_dendrogram_three_metrics_vertical.pdf\n")
cat("3. Dendrogram_sample_metadata.csv\n")

cat("\nAnalysis complete.\n")
