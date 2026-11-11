
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  COMPARE BURN-CONTROL DISTANCES AGAINST CONTROL_1-CONTROL_2 DISTANCE
#  All four metrics; years pooled
#This code uses the rZOTU object generated in other code within the GitHub folder
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(vegan)
library(phyloseq)
library(dplyr)

############################
# SETTINGS
############################

distance_metrics <- c(
  "jaccard",
  "bray",
  "unifrac",
  "wunifrac"
)

n_perm <- 9999

set.seed(123)

output_dir <- "Composition_Analysis_Results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

distance_names <- c(
  "jaccard" = "Jaccard",
  "bray" = "Bray-Curtis",
  "unifrac" = "Unweighted UniFrac",
  "wunifrac" = "Weighted UniFrac"
)

############################
# PREPARE COMMUNITY MATRIX
############################

comm_test <- as(
  otu_table(rZOTU),
  "matrix"
)

# Samples must be rows and taxa columns
if (taxa_are_rows(rZOTU)) {
  comm_test <- t(comm_test)
}

standardize_ids <- function(x) {
  tolower(gsub("-", ".", x))
}

rownames(comm_test) <- standardize_ids(
  rownames(comm_test)
)

############################
# PREPARE SITE LABELS
############################

meta_test <- as(
  sample_data(rZOTU),
  "data.frame"
)

rownames(meta_test) <- standardize_ids(
  rownames(meta_test)
)

# Match metadata to community matrix
meta_test <- meta_test[
  rownames(comm_test),
  ,
  drop = FALSE
]

if (!identical(
  rownames(comm_test),
  rownames(meta_test)
)) {
  stop("Community matrix and metadata are not aligned.")
}

site <- as.character(meta_test$Site)
names(site) <- rownames(meta_test)

if (anyNA(site)) {
  stop("Missing site labels.")
}

required_sites <- c(
  "Burn",
  "Control_1",
  "Control_2"
)

if (!all(required_sites %in% site)) {
  stop("One or more required sites are missing.")
}

cat("\nSample sizes by site:\n")
print(table(site))

############################
# FUNCTION: CALCULATE DISTANCE
############################

calculate_test_distance <- function(metric) {
  
  if (metric == "jaccard") {
    
    d <- vegan::vegdist(
      comm_test > 0,
      method = "jaccard",
      binary = TRUE
    )
    
  } else if (metric == "bray") {
    
    d <- vegan::vegdist(
      comm_test,
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
  
  # Align distance matrix with community matrix
  D <- as.matrix(d)
  
  rownames(D) <- standardize_ids(rownames(D))
  colnames(D) <- standardize_ids(colnames(D))
  
  sample_ids <- rownames(comm_test)
  
  if (
    !all(sample_ids %in% rownames(D)) ||
    !all(sample_ids %in% colnames(D))
  ) {
    stop("Distance matrix sample mismatch for ", metric)
  }
  
  D <- D[
    sample_ids,
    sample_ids,
    drop = FALSE
  ]
  
  if (any(!is.finite(D))) {
    stop("Non-finite distances detected for ", metric)
  }
  
  D
}

############################
# FUNCTION: MEAN BETWEEN-SITE DISTANCE
############################

mean_between <- function(group1, group2, labels, D) {
  
  i <- which(labels == group1)
  j <- which(labels == group2)
  
  if (length(i) == 0 || length(j) == 0) {
    stop(
      "No samples for comparison: ",
      group1,
      " versus ",
      group2
    )
  }
  
  mean(D[i, j, drop = FALSE])
}

############################
# FUNCTION: ONE-SIDED PERMUTATION TEST
############################

run_distance_comparison <- function(
    D,
    original_site,
    metric,
    comparison,
    n_perm = 9999,
    seed = 123) {
  
  # Comparison 1:
  # Burn-Control_1 < Control_1-Control_2
  #
  # Hold Control_1 fixed and permute Burn/Control_2 labels.
  
  if (comparison == "Burn-Control_1") {
    
    target_group1 <- "Burn"
    target_group2 <- "Control_1"
    
    reference_group1 <- "Control_1"
    reference_group2 <- "Control_2"
    
    permuted_groups <- c("Burn", "Control_2")
    
    # Comparison 2:
    # Burn-Control_2 < Control_1-Control_2
    #
    # Hold Control_2 fixed and permute Burn/Control_1 labels.
    
  } else if (comparison == "Burn-Control_2") {
    
    target_group1 <- "Burn"
    target_group2 <- "Control_2"
    
    reference_group1 <- "Control_1"
    reference_group2 <- "Control_2"
    
    permuted_groups <- c("Burn", "Control_1")
    
  } else {
    
    stop("Unknown comparison.")
  }
  
  ############################
  # OBSERVED DISTANCES
  ############################
  
  target_distance <- mean_between(
    target_group1,
    target_group2,
    original_site,
    D
  )
  
  reference_distance <- mean_between(
    reference_group1,
    reference_group2,
    original_site,
    D
  )
  
  observed_difference <- (
    target_distance - reference_distance
  )
  
  ############################
  # PERMUTATION TEST
  ############################
  
  # Permute labels only among the two exchangeable groups.
  # Keep the reference site's labels fixed.
  
  permuted_indices <- which(
    original_site %in% permuted_groups
  )
  
  perm_differences <- numeric(n_perm)
  
  set.seed(seed)
  
  for (b in seq_len(n_perm)) {
    
    perm_site <- original_site
    
    perm_site[permuted_indices] <- sample(
      original_site[permuted_indices],
      replace = FALSE
    )
    
    perm_target <- mean_between(
      target_group1,
      target_group2,
      perm_site,
      D
    )
    
    perm_reference <- mean_between(
      reference_group1,
      reference_group2,
      perm_site,
      D
    )
    
    perm_differences[b] <- (
      perm_target - perm_reference
    )
  }
  
  ############################
  # ONE-SIDED P-VALUE
  ############################
  
  # Alternative hypothesis:
  # target distance < reference distance
  
  p_value <- (
    sum(perm_differences <= observed_difference) + 1
  ) / (n_perm + 1)
  
  ############################
  # RETURN RESULTS
  ############################
  
  data.frame(
    Distance_Metric = unname(distance_names[metric]),
    Metric_Code = metric,
    Comparison = comparison,
    Target_Mean_Distance = target_distance,
    Reference_Mean_Distance = reference_distance,
    Difference_Target_Minus_Reference =
      observed_difference,
    P_value = p_value,
    N_Permutations = n_perm,
    stringsAsFactors = FALSE
  )
}

############################
# RUN BOTH TESTS FOR ALL FOUR METRICS
############################

all_results <- list()
result_index <- 1

for (metric in distance_metrics) {
  
  cat("\n========================================\n")
  cat("Distance metric:", distance_names[metric], "\n")
  cat("========================================\n")
  
  D <- calculate_test_distance(metric)
  
  # Ensure distance matrix and site labels are aligned
  site_metric <- site[rownames(D)]
  
  if (anyNA(site_metric)) {
    stop("Could not align site labels for ", metric)
  }
  
  # Test 1:
  # Is Burn-Control_1 distance smaller than Control_1-Control_2?
  
  result1 <- run_distance_comparison(
    D = D,
    original_site = site_metric,
    metric = metric,
    comparison = "Burn-Control_1",
    n_perm = n_perm,
    seed = 1000 + match(metric, distance_metrics) * 10 + 1
  )
  
  all_results[[result_index]] <- result1
  result_index <- result_index + 1
  
  # Test 2:
  # Is Burn-Control_2 distance smaller than Control_1-Control_2?
  
  result2 <- run_distance_comparison(
    D = D,
    original_site = site_metric,
    metric = metric,
    comparison = "Burn-Control_2",
    n_perm = n_perm,
    seed = 1000 + match(metric, distance_metrics) * 10 + 2
  )
  
  all_results[[result_index]] <- result2
  result_index <- result_index + 1
  
  cat("\nBurn-Control_1 comparison:\n")
  print(result1)
  
  cat("\nBurn-Control_2 comparison:\n")
  print(result2)
}

############################
# COMBINE RESULTS
############################

distance_comparison_table <- bind_rows(all_results)

############################
# BH CORRECTION ACROSS ALL EIGHT TESTS
############################

distance_comparison_table$P_adjusted_BH <- p.adjust(
  distance_comparison_table$P_value,
  method = "BH"
)

distance_comparison_table$Significant_BH_0.05 <- ifelse(
  distance_comparison_table$P_adjusted_BH < 0.05,
  "Yes",
  "No"
)

############################
# PRINT FINAL RESULTS
############################

cat("\n========================================\n")
cat("ALL DISTANCE COMPARISON RESULTS\n")
cat("========================================\n")

print(
  as.data.frame(distance_comparison_table),
  row.names = FALSE
)

############################
# SAVE RESULTS
############################

output_file <- file.path(
  output_dir,
  "Burn_Control_Distance_Comparisons_All_Metrics.csv"
)

write.csv(
  as.data.frame(distance_comparison_table),
  file = output_file,
  row.names = FALSE
)

cat("\nResults saved to:\n")
cat(output_file, "\n")

############################
# SUMMARY OF SIGNIFICANT RESULTS
############################

cat("\n========================================\n")
cat("SIGNIFICANT RESULTS AFTER BH CORRECTION\n")
cat("========================================\n")

significant_results <- distance_comparison_table %>%
  filter(
    P_adjusted_BH < 0.05,
    Difference_Target_Minus_Reference < 0
  )

if (nrow(significant_results) == 0) {
  
  cat("No comparisons were significant after BH correction.\n")
  
} else {
  
  print(
    as.data.frame(significant_results),
    row.names = FALSE
  )
}
