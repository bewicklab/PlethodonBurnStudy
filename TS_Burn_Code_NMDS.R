#===============================================================================
# MICROBIAL COMMUNITY COMPOSITION ANALYSES
# Four distance metrics: Jaccard, Bray-Curtis, Unweighted UniFrac, Weighted UniFrac
#
# Required object:
#   rZOTU - processed phyloseq object containing OTU table, taxonomy,
#           phylogenetic tree, and sample metadata. This object is generated in other codes within the GitHub folder.
#
# Analyses:
#   1. Overall PERMANOVA for each distance metric
#   2. Overall multivariate dispersion tests
#   3. NMDS plots for each distance metric and combined figure
#   4. Pairwise PERMANOVA among the three sites, with BH correction
#   5. Pairwise multivariate dispersion tests, with BH correction
#   6. NMDS stress table
#===============================================================================


#-------------------------------------------------------------------------------
# 1. PACKAGES AND SETTINGS
#-------------------------------------------------------------------------------

library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(patchwork)

distance_metrics <- c("jaccard", "bray", "unifrac", "wunifrac")

permutation_no <- 10000

set.seed(123)

output_dir <- "Composition_Analysis_Results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

site_levels <- c("Control_1", "Control_2", "Burn")

site_pairs <- list(
  c("Control_1", "Control_2"),
  c("Control_1", "Burn"),
  c("Control_2", "Burn")
)

cols <- c(
  "Control_1" = "#228A22",
  "Control_2" = "#00FE00",
  "Burn"      = "red"
)

year_shapes <- c(
  "2024" = 16,
  "2025" = 17
)

distance_names <- c(
  "jaccard" = "Jaccard",
  "bray" = "Bray-Curtis",
  "unifrac" = "Unweighted UniFrac",
  "wunifrac" = "Weighted UniFrac"
)


# Standardize sample IDs consistently throughout the analysis.

standardize_ids <- function(x) {
  tolower(gsub("-", ".", x))
}


#-------------------------------------------------------------------------------
# 2. HELPER FUNCTION: EXTRACT PERMANOVA RESULTS
#-------------------------------------------------------------------------------

# adonis2() can label the model row differently depending on the model.
# This function identifies the appropriate row and extracts the statistics.

extract_permanova <- function(fit) {
  
  fit_df <- as.data.frame(fit)
  
  possible_rows <- c("Model", "Site")
  
  model_row <- possible_rows[
    possible_rows %in% rownames(fit_df)
  ][1]
  
  if (is.na(model_row)) {
    stop(
      "Could not identify the model row in adonis2 output. Rows returned: ",
      paste(rownames(fit_df), collapse = ", ")
    )
  }
  
  required_columns <- c(
    "Df",
    "SumOfSqs",
    "R2",
    "F",
    "Pr(>F)"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(fit_df)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing expected PERMANOVA columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  data.frame(
    Df = as.numeric(fit_df[model_row, "Df"]),
    SumOfSqs = as.numeric(fit_df[model_row, "SumOfSqs"]),
    R2 = as.numeric(fit_df[model_row, "R2"]),
    F_statistic = as.numeric(fit_df[model_row, "F"]),
    P_value = as.numeric(fit_df[model_row, "Pr(>F)"])
  )
}


#-------------------------------------------------------------------------------
# 3. PREPARE COMMUNITY MATRIX AND METADATA
#-------------------------------------------------------------------------------

if (!exists("rZOTU")) {
  stop("The phyloseq object 'rZOTU' was not found in the R environment.")
}

comm <- as(
  otu_table(rZOTU),
  "matrix"
)

# vegan expects samples in rows and taxa in columns.

if (taxa_are_rows(rZOTU)) {
  comm <- t(comm)
}

rownames(comm) <- standardize_ids(
  rownames(comm)
)

if (anyDuplicated(rownames(comm))) {
  stop(
    "Sample IDs are duplicated after standardization in the community matrix."
  )
}

meta <- as(
  sample_data(rZOTU),
  "data.frame"
)

rownames(meta) <- standardize_ids(
  rownames(meta)
)

if (anyDuplicated(rownames(meta))) {
  stop(
    "Sample IDs are duplicated after standardization in the metadata."
  )
}

if (!all(rownames(comm) %in% rownames(meta))) {
  stop("Some community samples are missing from the metadata.")
}

# Reorder metadata to match the community matrix exactly.

meta <- meta[
  rownames(comm),
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(comm), rownames(meta))
)


#-------------------------------------------------------------------------------
# 4. DEFINE SITE AND YEAR GROUPS
#-------------------------------------------------------------------------------

if (!"Site" %in% colnames(meta)) {
  stop("The metadata does not contain a 'Site' column.")
}

meta$Site <- factor(
  as.character(meta$Site),
  levels = site_levels
)

if (anyNA(meta$Site)) {
  stop("Some samples have missing or unmatched Site values.")
}

if ("Year" %in% colnames(meta)) {
  
  year_values <- as.character(meta$Year)
  
} else if ("Season" %in% colnames(meta)) {
  
  year_values <- sub(
    ".*_(20[0-9]{2})$",
    "\\1",
    as.character(meta$Season)
  )
  
} else {
  
  stop("Could not find a Year or Season column in the metadata.")
  
}

year_groups <- factor(
  year_values,
  levels = c("2024", "2025")
)

if (anyNA(year_groups)) {
  stop("Some samples do not have a valid 2024 or 2025 year label.")
}


#-------------------------------------------------------------------------------
# 5. CALCULATE AND ALIGN ALL FOUR DISTANCE MATRICES ONCE
#-------------------------------------------------------------------------------

calculate_distance <- function(metric) {
  
  if (metric == "jaccard") {
    
    d <- vegan::vegdist(
      comm > 0,
      method = "jaccard",
      binary = TRUE
    )
    
  } else if (metric == "bray") {
    
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
  
  if (is.null(rownames(dmat)) || is.null(colnames(dmat))) {
    stop("Distance matrix has no sample names for metric: ", metric)
  }
  
  rownames(dmat) <- standardize_ids(
    rownames(dmat)
  )
  
  colnames(dmat) <- standardize_ids(
    colnames(dmat)
  )
  
  if (
    anyDuplicated(rownames(dmat)) ||
    anyDuplicated(colnames(dmat))
  ) {
    stop(
      "Distance matrix sample IDs are duplicated after standardization: ",
      metric
    )
  }
  
  sample_ids <- rownames(comm)
  
  if (
    !all(sample_ids %in% rownames(dmat)) ||
    !all(sample_ids %in% colnames(dmat))
  ) {
    stop(
      "Distance matrix samples do not match the community matrix: ",
      metric
    )
  }
  
  # Explicitly align both dimensions to the community matrix.
  
  dmat <- dmat[
    sample_ids,
    sample_ids,
    drop = FALSE
  ]
  
  stopifnot(
    identical(rownames(dmat), sample_ids),
    identical(colnames(dmat), sample_ids)
  )
  
  as.dist(dmat)
}

distance_matrices <- setNames(
  lapply(distance_metrics, calculate_distance),
  distance_metrics
)


#-------------------------------------------------------------------------------
# 6. NMDS PLOT FUNCTION
#-------------------------------------------------------------------------------

make_nmds_plot <- function(metric, nmds_result, scores_nmds) {
  
  plot_data <- data.frame(
    Sample = rownames(scores_nmds),
    NMDS1 = scores_nmds[, 1],
    NMDS2 = scores_nmds[, 2],
    Site = meta$Site,
    Year = year_groups,
    stringsAsFactors = FALSE
  )
  
  plot_data$Site <- factor(
    plot_data$Site,
    levels = site_levels
  )
  
  plot_data$Year <- factor(
    plot_data$Year,
    levels = c("2024", "2025")
  )
  
  # Calculate one convex hull per site, when at least
  # three distinct points exist.
  
  hull_list <- lapply(site_levels, function(site) {
    
    pts <- plot_data[
      plot_data$Site == site,
      c("NMDS1", "NMDS2"),
      drop = FALSE
    ]
    
    if (
      nrow(pts) < 3 ||
      nrow(unique(pts)) < 3
    ) {
      return(NULL)
    }
    
    hull_indices <- chull(
      pts$NMDS1,
      pts$NMDS2
    )
    
    hull_indices <- c(
      hull_indices,
      hull_indices[1]
    )
    
    data.frame(
      NMDS1 = pts$NMDS1[hull_indices],
      NMDS2 = pts$NMDS2[hull_indices],
      Site = site
    )
  })
  
  hull_data <- bind_rows(hull_list)
  
  hull_data$Site <- factor(
    hull_data$Site,
    levels = site_levels
  )
  
  ggplot() +
    
    geom_polygon(
      data = hull_data,
      aes(
        x = NMDS1,
        y = NMDS2,
        group = Site,
        fill = Site,
        color = Site
      ),
      alpha = 0.20,
      linewidth = 0.8
    ) +
    
    geom_point(
      data = plot_data,
      aes(
        x = NMDS1,
        y = NMDS2,
        color = Site,
        shape = Year
      ),
      size = 2.8
    ) +
    
    scale_color_manual(
      values = cols,
      breaks = names(cols),
      labels = c("Control 1", "Control 2", "Burn"),
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = cols,
      breaks = names(cols),
      labels = c("Control 1", "Control 2", "Burn"),
      drop = FALSE
    ) +
    
    scale_shape_manual(
      values = year_shapes,
      breaks = c("2024", "2025"),
      drop = FALSE
    ) +
    
    labs(
      title = unname(distance_names[metric]),
      subtitle = sprintf(
        "Stress = %.4f",
        nmds_result$stress
      ),
      x = "NMDS1",
      y = "NMDS2",
      color = "Site",
      fill = "Site",
      shape = "Year"
    ) +
    
    theme_classic(base_size = 12) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5
      ),
      axis.title = element_text(size = 11),
      legend.position = "right"
    )
}


#-------------------------------------------------------------------------------
# 7. OVERALL PERMANOVA, OVERALL DISPERSION, AND NMDS
#-------------------------------------------------------------------------------

nmds_plots <- list()
nmds_results <- list()
permanova_results <- list()
dispersion_results <- list()

for (metric in distance_metrics) {
  
  cat("\n========================================\n")
  cat("Distance metric:", metric, "\n")
  cat("========================================\n")
  
  d <- distance_matrices[[metric]]
  
  # Overall PERMANOVA: do sites differ in community composition?
  
  permanova_fit <- vegan::adonis2(
    d ~ Site,
    data = meta,
    permutations = permutation_no
  )
  
  cat("\nOverall PERMANOVA:\n")
  print(permanova_fit)
  
  # FIXED: robust extraction of the PERMANOVA model row.
  
  permanova_stats <- extract_permanova(
    permanova_fit
  )
  
  permanova_results[[metric]] <- data.frame(
    Distance_Metric = unname(distance_names[metric]),
    Metric_Code = metric,
    Term = "Site",
    permanova_stats
  )
  
  # Overall homogeneity of multivariate dispersion.
  
  bd <- vegan::betadisper(
    d,
    group = meta$Site
  )
  
  bd_anova <- anova(bd)
  
  cat("\nOverall betadisper:\n")
  print(bd_anova)
  
  dispersion_results[[metric]] <- data.frame(
    Distance_Metric = unname(distance_names[metric]),
    Metric_Code = metric,
    Df = bd_anova[1, "Df"],
    F_statistic = bd_anova[1, "F value"],
    P_value = bd_anova[1, "Pr(>F)"]
  )
  
  # Two-dimensional NMDS.
  
  set.seed(123)
  
  nmds <- vegan::metaMDS(
    d,
    k = 2,
    trymax = 500,
    autotransform = FALSE,
    trace = FALSE
  )
  
  cat(
    "\nNMDS stress for ",
    unname(distance_names[metric]),
    ": ",
    nmds$stress,
    "\n",
    sep = ""
  )
  
  scores_nmds <- as.matrix(
    vegan::scores(
      nmds,
      display = "sites",
      choices = 1:2
    )
  )
  
  if (ncol(scores_nmds) < 2) {
    stop(
      "NMDS returned fewer than two dimensions for ",
      metric
    )
  }
  
  if (!all(rownames(comm) %in% rownames(scores_nmds))) {
    stop(
      "NMDS scores could not be matched to samples for ",
      metric
    )
  }
  
  scores_nmds <- scores_nmds[
    rownames(comm),
    1:2,
    drop = FALSE
  ]
  
  nmds_plot <- make_nmds_plot(
    metric,
    nmds,
    scores_nmds
  )
  
  nmds_plots[[metric]] <- nmds_plot
  nmds_results[[metric]] <- nmds
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0("NMDS_", metric, ".pdf")
    ),
    plot = nmds_plot,
    width = 8,
    height = 6,
    units = "in"
  )
}


#-------------------------------------------------------------------------------
# 8. SAVE OVERALL PERMANOVA AND DISPERSION RESULTS
#-------------------------------------------------------------------------------

permanova_table <- bind_rows(
  permanova_results
)

dispersion_table <- bind_rows(
  dispersion_results
)

cat("\nOVERALL PERMANOVA RESULTS\n")
print(permanova_table)

write.csv(
  permanova_table,
  file.path(
    output_dir,
    "PERMANOVA_all_distance_metrics.csv"
  ),
  row.names = FALSE
)

cat("\nOVERALL BETADISPER RESULTS\n")
print(dispersion_table)

write.csv(
  dispersion_table,
  file.path(
    output_dir,
    "Betadisper_all_distance_metrics.csv"
  ),
  row.names = FALSE
)


#-------------------------------------------------------------------------------
# 9. COMBINE AND SAVE FOUR-PANEL NMDS FIGURE
#-------------------------------------------------------------------------------

combined_nmds <- (
  nmds_plots[["jaccard"]] |
    nmds_plots[["bray"]]
) / (
  nmds_plots[["unifrac"]] |
    nmds_plots[["wunifrac"]]
)

combined_nmds <- combined_nmds +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_nmds)

ggsave(
  filename = file.path(
    output_dir,
    "NMDS_all_four_distance_metrics.pdf"
  ),
  plot = combined_nmds,
  width = 14,
  height = 11,
  units = "in"
)


#-------------------------------------------------------------------------------
# 10. SAVE NMDS STRESS VALUES
#-------------------------------------------------------------------------------

stress_table <- data.frame(
  Distance_Metric = unname(
    distance_names[names(nmds_results)]
  ),
  Metric_Code = names(nmds_results),
  NMDS_Stress = vapply(
    nmds_results,
    function(x) x$stress,
    numeric(1)
  )
)

cat("\nNMDS STRESS VALUES\n")
print(stress_table)

write.csv(
  stress_table,
  file.path(
    output_dir,
    "NMDS_stress_all_distance_metrics.csv"
  ),
  row.names = FALSE
)


#-------------------------------------------------------------------------------
# 11. PAIRWISE PERMANOVA
#-------------------------------------------------------------------------------

pairwise_permanova_results <- list()

result_index <- 1

for (metric in distance_metrics) {
  
  dmat <- as.matrix(
    distance_matrices[[metric]]
  )
  
  for (i in seq_along(site_pairs)) {
    
    pair <- site_pairs[[i]]
    
    keep <- as.character(meta$Site) %in% pair
    
    meta_sub <- meta[
      keep,
      ,
      drop = FALSE
    ]
    
    meta_sub$Site <- factor(
      as.character(meta_sub$Site),
      levels = pair
    )
    
    sample_ids <- rownames(meta_sub)
    
    d_sub <- as.dist(
      dmat[
        sample_ids,
        sample_ids,
        drop = FALSE
      ]
    )
    
    # Use a reproducible, distinct seed for each comparison.
    
    set.seed(
      123 + match(metric, distance_metrics) * 100 + i
    )
    
    pair_fit <- vegan::adonis2(
      d_sub ~ Site,
      data = meta_sub,
      permutations = permutation_no
    )
    
    cat(
      "\nPairwise PERMANOVA:",
      metric,
      "-",
      pair[1],
      "vs",
      pair[2],
      "\n"
    )
    
    print(pair_fit)
    
    # FIXED: use the same robust extraction function.
    
    pair_stats <- extract_permanova(
      pair_fit
    )
    
    pairwise_permanova_results[[result_index]] <- data.frame(
      Distance_Metric = unname(distance_names[metric]),
      Metric_Code = metric,
      Site_1 = pair[1],
      Site_2 = pair[2],
      Comparison = paste(
        pair[1],
        "vs",
        pair[2]
      ),
      N_Site_1 = sum(meta_sub$Site == pair[1]),
      N_Site_2 = sum(meta_sub$Site == pair[2]),
      N_Samples = nrow(meta_sub),
      pair_stats
    )
    
    result_index <- result_index + 1
  }
}

pairwise_permanova_table <- bind_rows(
  pairwise_permanova_results
) %>%
  group_by(Metric_Code) %>%
  mutate(
    P_adjusted_BH = p.adjust(
      P_value,
      method = "BH"
    ),
    Significant_BH_0.05 = ifelse(
      P_adjusted_BH < 0.05,
      "Yes",
      "No"
    )
  ) %>%
  ungroup()

cat("\nALL PAIRWISE PERMANOVA RESULTS\n")

print(
  as.data.frame(pairwise_permanova_table),
  row.names = FALSE
)

write.csv(
  as.data.frame(pairwise_permanova_table),
  file.path(
    output_dir,
    "Pairwise_PERMANOVA_all_distance_metrics.csv"
  ),
  row.names = FALSE
)


#-------------------------------------------------------------------------------
# 12. PAIRWISE MULTIVARIATE DISPERSION TESTS
#-------------------------------------------------------------------------------

pairwise_dispersion_results <- list()

result_index <- 1

for (metric in distance_metrics) {
  
  dmat <- as.matrix(
    distance_matrices[[metric]]
  )
  
  for (i in seq_along(site_pairs)) {
    
    pair <- site_pairs[[i]]
    
    keep <- as.character(meta$Site) %in% pair
    
    meta_sub <- meta[
      keep,
      ,
      drop = FALSE
    ]
    
    sample_ids <- rownames(meta_sub)
    
    d_sub <- as.dist(
      dmat[
        sample_ids,
        sample_ids,
        drop = FALSE
      ]
    )
    
    group_sub <- droplevels(
      factor(
        as.character(meta_sub$Site),
        levels = pair
      )
    )
    
    # Distances of samples from their group centroids.
    
    bd <- vegan::betadisper(
      d_sub,
      group = group_sub
    )
    
    set.seed(
      123 + match(metric, distance_metrics) * 100 + i
    )
    
    dispersion_test <- vegan::permutest(
      bd,
      permutations = permutation_no
    )
    
    test_table <- dispersion_test$tab
    
    mean_distances <- tapply(
      bd$distances,
      group_sub,
      mean,
      na.rm = TRUE
    )
    
    pairwise_dispersion_results[[result_index]] <- data.frame(
      Distance_Metric = unname(distance_names[metric]),
      Metric_Code = metric,
      Site_1 = pair[1],
      Site_2 = pair[2],
      Comparison = paste(
        pair[1],
        "vs",
        pair[2]
      ),
      N_Site_1 = sum(group_sub == pair[1]),
      N_Site_2 = sum(group_sub == pair[2]),
      N_Samples = length(group_sub),
      Mean_Distance_Site_1 = unname(
        mean_distances[pair[1]]
      ),
      Mean_Distance_Site_2 = unname(
        mean_distances[pair[2]]
      ),
      F_statistic = as.numeric(
        test_table[1, "F"]
      ),
      P_value = as.numeric(
        test_table[1, "Pr(>F)"]
      )
    )
    
    cat(
      "\nPairwise dispersion:",
      metric,
      "-",
      pair[1],
      "vs",
      pair[2],
      "\n"
    )
    
    print(
      pairwise_dispersion_results[[result_index]]
    )
    
    result_index <- result_index + 1
  }
}

pairwise_dispersion_table <- bind_rows(
  pairwise_dispersion_results
) %>%
  group_by(Metric_Code) %>%
  mutate(
    P_adjusted_BH = p.adjust(
      P_value,
      method = "BH"
    ),
    Significant_BH_0.05 = ifelse(
      P_adjusted_BH < 0.05,
      "Yes",
      "No"
    )
  ) %>%
  ungroup()

cat("\nALL PAIRWISE DISPERSION RESULTS\n")

print(
  as.data.frame(pairwise_dispersion_table),
  row.names = FALSE
)

write.csv(
  as.data.frame(pairwise_dispersion_table),
  file.path(
    output_dir,
    "Pairwise_Betadisper_all_distance_metrics.csv"
  ),
  row.names = FALSE
)


#-------------------------------------------------------------------------------
# 13. FINAL SUMMARY
#-------------------------------------------------------------------------------

cat("\n========================================\n")
cat("ALL FOUR DISTANCE METRICS COMPLETED\n")
cat("========================================\n")

cat(
  "Output directory:",
  output_dir,
  "\n\n"
)

cat("Files saved:\n")

cat("1. NMDS_all_four_distance_metrics.pdf\n")
cat("2. NMDS_jaccard.pdf\n")
cat("3. NMDS_bray.pdf\n")
cat("4. NMDS_unifrac.pdf\n")
cat("5. NMDS_wunifrac.pdf\n")
cat("6. PERMANOVA_all_distance_metrics.csv\n")
cat("7. Betadisper_all_distance_metrics.csv\n")
cat("8. NMDS_stress_all_distance_metrics.csv\n")
cat("9. Pairwise_PERMANOVA_all_distance_metrics.csv\n")
cat("10. Pairwise_Betadisper_all_distance_metrics.csv\n")

cat("\nAnalysis complete.\n")