############################################################
# STREAMLINED MICROBIAL TAXON + FAITH'S PD ACCUMULATION
# Outputs retained from the original script.
############################################################

# 1. PACKAGES AND SETTINGS -----------------------------------
required <- c("qiime2R", "phyloseq", "dplyr", "ggplot2", "patchwork", "ape")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "))
invisible(lapply(required, library, character.only = TRUE))

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

set.seed(123)
n_boot <- 999L
n_perm <- 999L
n_subsample <- 20L
alpha <- 0.05
site_levels <- c("Control_1", "Control_2", "Burn")
site_colors <- c(Control_1 = "#228A22", Control_2 = "#00FE00", Burn = "red")
site_labels <- c(Control_1 = "Control 1", Control_2 = "Control 2", Burn = "Burn")
site_pairs <- combn(site_levels, 2, simplify = FALSE)
panel_definitions <- list(
  "2024" = "2024",
  "2025" = "2025",
  "2024 + 2025" = c("2024", "2025")
)

# 2. IMPORT AND PREPROCESS QIIME DATA -------------------------
# Run from the directory containing table.qza, taxonomy.qza,
# tree.qza, and Coweeta_TEWL_Seasons_Sites.csv.
meta <- read.csv("Coweeta_TEWL_Seasons_Sites.csv", stringsAsFactors = FALSE)
meta <- meta %>%
  filter(!is.na(QIIME_ID_aligned), QIIME_ID_aligned != "") %>%
  distinct(QIIME_ID_aligned, .keep_all = TRUE)
rownames(meta) <- tolower(gsub("-", ".", meta$QIIME_ID_aligned))

ps <- qiime2R::qza_to_phyloseq(
  features = "table.qza",
  taxonomy = "taxonomy.qza",
  tree = "tree.qza"
)
sample_names(ps) <- tolower(gsub("-", ".", sample_names(ps)))
common_ids <- intersect(sample_names(ps), rownames(meta))
if (!length(common_ids)) stop("No QIIME sample IDs match the metadata IDs.")
ps <- prune_samples(common_ids, ps)
meta <- meta[common_ids, , drop = FALSE]
sample_data(ps) <- sample_data(meta)
stopifnot(identical(sample_names(ps), rownames(data.frame(sample_data(ps)))))

# Remove rare taxa before blank subtraction, as in the original.
ps <- prune_taxa(taxa_sums(ps) > 10, ps)
blank_ids <- sample_names(ps)[as.character(sample_data(ps)$Site) == "Blank"]
if (!length(blank_ids)) stop("No Blank control samples detected.")
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
blank_max <- apply(otu[, blank_ids, drop = FALSE], 1, max, na.rm = TRUE)
blank_max[!is.finite(blank_max)] <- 0
otu <- sweep(otu, 1, blank_max, "-")
otu[otu < 0] <- 0
otu_table(ps) <- otu_table(otu, taxa_are_rows = TRUE)
ps <- prune_samples(as.character(sample_data(ps)$Site) != "Blank", ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Remove unwanted taxonomic groups.
tax <- as.data.frame(tax_table(ps))
keep_taxa <- !(as.character(tax$Kingdom) %in% "d__Eukaryota") &
  !(as.character(tax$Phylum) %in% c(
    "Chytridiomycota", "Ascomycota", "Chlorophyta", "Apicomplexa",
    "Ciliophora", "Phragmoplastophyta", "Zoopagomycota"
  )) & as.character(tax$Family) != "Chloroplast" &
  as.character(tax$Genus) != "Mitochondria"
keep_taxa[is.na(keep_taxa)] <- FALSE
ps <- prune_taxa(keep_taxa, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Exclude the same sites, samples, seasons, and missing-coordinate records.
sd <- data.frame(sample_data(ps))
exclude_sites <- c("Control", "2750", "2850", "EF", "PN", "Positive_control")
exclude_titles <- paste0("TS.31.", 5:11)
exclude_seasons <- c("Summer_2024", "Summer_2025", "Fall_2024")
keep_samples <- !as.character(sd$Site) %in% exclude_sites &
  !as.character(sd$title) %in% exclude_titles &
  !as.character(sd$Season) %in% exclude_seasons &
  !is.na(sd$lat_7_Collection_locatio)
keep_samples[is.na(keep_samples)] <- FALSE
ps <- prune_samples(keep_samples, ps)
if (!nsamples(ps)) stop("No samples remain after preprocessing.")

# Rarefy to 12,000 reads, without replacement.
rZOTU <- rarefy_even_depth(ps, sample.size = 12000, rngseed = 123, replace = FALSE, verbose = TRUE)

# Standardize site names and create Site x Year metadata.
sd <- data.frame(sample_data(rZOTU))
sd$Site <- as.character(sd$Site)
sd$Site[sd$Site == "31_burn"] <- "Burn"
sd$Site[sd$Site == "32_no_burn"] <- "Control_2"
sd$Site[sd$Site == "2950"] <- "Control_1"
sd$Year <- sub("^.*_", "", as.character(sd$Season))
sd$SitexYear <- paste(sd$Site, sd$Year, sep = "_")
sample_data(rZOTU) <- sample_data(sd)

# One shared sample-by-taxon presence/absence matrix and metadata table.
otu <- as(otu_table(rZOTU), "matrix")
if (taxa_are_rows(rZOTU)) otu <- t(otu)
metadata <- data.frame(
  Sample = rownames(otu),
  Site = as.character(sample_data(rZOTU)$Site),
  Season = as.character(sample_data(rZOTU)$Season),
  Year = sub("^.*_", "", as.character(sample_data(rZOTU)$Season)),
  stringsAsFactors = FALSE
)
keep <- metadata$Site %in% site_levels & metadata$Year %in% c("2024", "2025")
metadata <- metadata[keep, , drop = FALSE]
otu <- otu[metadata$Sample, , drop = FALSE]
community <- 1L * (otu > 0)
community <- community[, colSums(community) > 0, drop = FALSE]
if (!nrow(community) || !ncol(community)) stop("The retained community matrix is empty.")

# 3. SHARED PLOTTING HELPERS ---------------------------------
make_theme <- function() {
  theme_classic(base_size = 16) + theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 21, face = "bold"),
    legend.text = element_text(size = 19),
    legend.key.size = grid::unit(1.3, "cm"),
    legend.key.width = grid::unit(1.3, "cm"),
    legend.spacing.y = grid::unit(0.4, "cm"),
    legend.position = "right"
  )
}

combine_panels <- function(results) {
  (patchwork::wrap_plots(
    A = results[[1]]$plot, B = results[[2]]$plot, C = results[[3]]$plot,
    design = "AB\nC#", guides = "collect", widths = c(1, 1), heights = c(1, 1)
  ) & theme(legend.position = "right"))
}

# 4. NON-PHYLOGENETIC TAXON ACCUMULATION ----------------------
# Expected richness after k samples, using the hypergeometric
# probability that each taxon is absent from the selected samples.
expected_richness <- function(x, k) {
  n <- nrow(x)
  freq <- colSums(x > 0)
  absent_n <- n - freq
  p_absent <- numeric(length(absent_n))
  possible <- absent_n >= k
  p_absent[possible] <- exp(lchoose(absent_n[possible], k) - lchoose(n, k))
  sum(1 - p_absent)
}

bootstrap_curve <- function(x, effort, metric = c("richness", "pd")) {
  metric <- match.arg(metric)
  n <- nrow(x)
  out <- matrix(NA_real_, nrow = n_boot, ncol = effort)
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    xb <- x[idx, , drop = FALSE]
    xb <- xb[sample(seq_len(n)), , drop = FALSE]
    if (metric == "richness") {
      out[b, ] <- vapply(seq_len(effort), function(k) expected_richness(xb, k), numeric(1))
    } else {
      seen_taxa <- logical(ncol(xb))
      seen_edges <- logical(pd_n_edges)
      value <- 0
      for (k in seq_len(effort)) {
        new_taxa <- which(xb[k, ] > 0 & !seen_taxa)
        if (length(new_taxa)) {
          seen_taxa[new_taxa] <- TRUE
          for (j in new_taxa) {
            node <- pd_tip_ids[j]
            repeat {
              if (node < 1L || node > pd_max_node) break
              edge_id <- pd_parent_edge[node]
              if (edge_id == 0L || seen_edges[edge_id]) break
              seen_edges[edge_id] <- TRUE
              value <- value + pd_edge_lengths[edge_id]
              node <- pd_parent_node[node]
            }
          }
        }
        out[b, k] <- value
      }
    }
  }
  out
}

make_accumulation_panel <- function(years, title, metric = c("richness", "pd")) {
  metric <- match.arg(metric)
  md <- metadata[metadata$Year %in% years, , drop = FALSE]
  counts <- table(factor(md$Site, levels = site_levels))
  if (any(counts == 0)) stop("At least one site has no samples in panel: ", title)
  effort <- min(counts)
  cat("\n", title, " sample counts: ", paste(names(counts), counts, sep = "=", collapse = "; "),
      " | common effort: ", effort, "\n", sep = "")
  summaries <- vector("list", length(site_levels))
  for (i in seq_along(site_levels)) {
    site <- site_levels[i]
    x <- community[md$Sample[md$Site == site], , drop = FALSE]
    if (metric == "pd") x <- x[, pd_taxa, drop = FALSE]
    boot <- bootstrap_curve(x, effort, metric)
    summaries[[i]] <- data.frame(
      Panel = title, Site = site, Samples = seq_len(effort),
      Mean = colMeans(boot),
      Lower_95_CI = apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE),
      Upper_95_CI = apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
    rm(boot); gc(FALSE)
  }
  summary <- do.call(rbind, summaries)
  names(summary)[names(summary) == "Mean"] <- if (metric == "richness") "Mean_Richness" else "Mean_PD"
  yvar <- if (metric == "richness") "Mean_Richness" else "Mean_PD"
  ylab <- if (metric == "richness") "Cumulative number of taxa" else "Faith's phylogenetic diversity"
  p <- ggplot(summary, aes(x = Samples, y = .data[[yvar]], color = Site, fill = Site, group = Site)) +
    geom_ribbon(aes(ymin = Lower_95_CI, ymax = Upper_95_CI), alpha = 0.20, color = NA) +
    geom_line(linewidth = 1) + geom_point(size = 1.5) +
    scale_color_manual(values = site_colors, labels = site_labels, drop = FALSE) +
    scale_fill_manual(values = site_colors, labels = site_labels, drop = FALSE) +
    scale_x_continuous(breaks = pretty(c(1, effort), n = 6), limits = c(1, effort)) +
    labs(title = title, x = "Number of samples", y = ylab, color = "Site", fill = "Site") +
    make_theme()
  list(plot = p, summary = summary)
}

run_curve_set <- function(metric) {
  results <- lapply(names(panel_definitions), function(nm) {
    make_accumulation_panel(panel_definitions[[nm]], nm, metric)
  })
  names(results) <- names(panel_definitions)
  fig <- combine_panels(results)
  print(fig)
  if (metric == "richness") {
    ggsave("Nonphylogenetic_Taxon_Accumulation_2024_2025_Both.pdf", fig, width = 14, height = 12, units = "in")
    ggsave("Nonphylogenetic_Taxon_Accumulation_2024_2025_Both.png", fig, width = 14, height = 12, units = "in", dpi = 300)
    write.csv(do.call(rbind, lapply(results, `[[`, "summary")),
              "Nonphylogenetic_Taxon_Accumulation_Summary_All_Years.csv", row.names = FALSE)
  } else {
    ggsave("Faith_PD_Accumulation_2024_2025_Both.pdf", fig, width = 14, height = 12, units = "in")
    ggsave("Faith_PD_Accumulation_2024_2025_Both.png", fig, width = 14, height = 12, units = "in", dpi = 300)
    write.csv(do.call(rbind, lapply(results, `[[`, "summary")),
              "Faith_PD_Accumulation_Summary_All_Years.csv", row.names = FALSE)
  }
  results
}

richness_curves <- run_curve_set("richness")

# 5. GENERIC PAIRWISE PERMUTATION TEST ------------------------
# Retains original design: 20 repeated equal-effort subsamples
# per observed/permuted statistic; 999 label permutations;
# absolute mean signed difference; Holm correction over 9 tests.
mean_metric_difference <- function(ids1, ids2, effort, sample_metric, n_sub = n_subsample) {
  diffs <- numeric(n_sub)
  for (i in seq_len(n_sub)) {
    a <- sample(ids1, effort, replace = FALSE)
    b <- sample(ids2, effort, replace = FALSE)
    diffs[i] <- mean(sample_metric(a)) - mean(sample_metric(b))
  }
  mean(diffs)
}

run_pairwise_tests <- function(metric_name, sample_metric, sample_ids_available) {
  all_results <- list()
  out_i <- 1L
  for (panel_name in names(panel_definitions)) {
    years <- panel_definitions[[panel_name]]
    md <- metadata[metadata$Year %in% years & metadata$Site %in% site_levels &
                     metadata$Sample %in% sample_ids_available, , drop = FALSE]
    md <- md[!duplicated(md$Sample), , drop = FALSE]
    ids_by_site <- setNames(lapply(site_levels, function(s) as.character(md$Sample[md$Site == s])), site_levels)
    counts <- lengths(ids_by_site)
    if (any(counts == 0)) stop("A site has no usable samples in panel: ", panel_name)
    effort <- min(counts)
    cat("\n", metric_name, " | ", panel_name, " | counts: ",
        paste(names(counts), counts, sep = "=", collapse = "; "),
        " | common effort: ", effort, "\n", sep = "")
    for (pair in site_pairs) {
      s1 <- pair[1]; s2 <- pair[2]
      ids1 <- ids_by_site[[s1]]; ids2 <- ids_by_site[[s2]]
      n1 <- length(ids1); n2 <- length(ids2)
      pooled <- c(ids1, ids2)
      labels <- c(rep(s1, n1), rep(s2, n2))
      observed_diff <- mean_metric_difference(ids1, ids2, effort, sample_metric)
      observed <- abs(observed_diff)
      null <- numeric(n_perm)
      progress <- max(1L, ceiling(n_perm / 10))
      start <- Sys.time()
      for (b in seq_len(n_perm)) {
        perm_labels <- sample(labels)
        p1 <- pooled[perm_labels == s1]
        p2 <- pooled[perm_labels == s2]
        null[b] <- abs(mean_metric_difference(p1, p2, effort, sample_metric))
        if (b %% progress == 0L || b == n_perm) {
          elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
          cat(sprintf("  %s vs %s: %d/%d permutations (%.0f%%); %.1f sec\n",
                      s1, s2, b, n_perm, 100 * b / n_perm, elapsed))
        }
      }
      p <- (1 + sum(null >= observed)) / (n_perm + 1)
      all_results[[out_i]] <- data.frame(
        Panel = panel_name, Site_1 = s1, Site_2 = s2,
        Available_N_1 = n1, Available_N_2 = n2,
        Common_Effort = effort, Subsamples_Per_Statistic = n_subsample,
        Mean_Difference = observed_diff, Test_Statistic = observed,
        Permutations = n_perm, P_Value = p, stringsAsFactors = FALSE
      )
      out_i <- out_i + 1L
    }
  }
  ans <- do.call(rbind, all_results)
  rownames(ans) <- NULL
  ans$Site_1_Label <- unname(site_labels[ans$Site_1])
  ans$Site_2_Label <- unname(site_labels[ans$Site_2])
  ans$P_Adjusted_Holm <- p.adjust(ans$P_Value, method = "holm")
  ans$Significant_Holm <- ans$P_Adjusted_Holm < alpha
  print(ans, digits = 8, row.names = FALSE)
  ans
}

# Richness permutation test: richness is the number of unique taxa
# in the selected samples.
taxa_by_sample <- lapply(seq_len(nrow(community)), function(i) which(community[i, ] > 0))
names(taxa_by_sample) <- rownames(community)
richness_sample_metric <- function(ids) length(unique(unlist(taxa_by_sample[ids], use.names = FALSE)))
richness_tests <- run_pairwise_tests("Non-phylogenetic richness", richness_sample_metric, names(taxa_by_sample))
significant_richness <- subset(richness_tests, Significant_Holm)
openxlsx::write.xlsx(richness_tests, "Nonphylogenetic_MaxCommonEffort_Pairwise_Tests.xlsx", overwrite = TRUE, rowNames = FALSE)
openxlsx::write.xlsx(significant_richness, "Nonphylogenetic_MaxCommonEffort_Significant_Tests.xlsx", overwrite = TRUE, rowNames = FALSE)

# 6. TREE DIAGNOSTICS AND FAITH'S PD SETUP --------------------
tree <- phy_tree(rZOTU, errorIfNULL = FALSE)
if (is.null(tree) || !inherits(tree, "phylo")) stop("rZOTU must contain an ape phylogenetic tree.")
if (is.null(tree$edge.length)) stop("Faith's PD requires branch lengths.")
if (anyNA(tree$edge.length) || any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
  stop("Tree branch lengths must be finite, non-missing, and nonnegative.")
}
if (anyDuplicated(tree$tip.label)) stop("The phylogenetic tree contains duplicated tip labels.")
cat("\nTree tips: ", length(tree$tip.label), " | branches: ", nrow(tree$edge), "\n", sep = "")
cat("Community taxa: ", ncol(community), "\n", sep = "")
cat("Matching taxa: ", length(intersect(colnames(community), tree$tip.label)), "\n", sep = "")
cat("Community taxa absent from tree (first 20):\n")
print(head(setdiff(colnames(community), tree$tip.label), 20))
cat("Tree tips absent from community (first 20):\n")
print(head(setdiff(tree$tip.label, colnames(community)), 20))
plot(tree, cex = 0.4, no.margin = TRUE)
ape::write.tree(tree, file = "rZOTU_phylogenetic_tree.nwk")

pd_taxa <- intersect(colnames(community), tree$tip.label)
if (!length(pd_taxa)) stop("No matching taxa between community matrix and tree.")
pd_community <- community[, pd_taxa, drop = FALSE]
pd_tree <- ape::keep.tip(tree, pd_taxa)
pd_community <- pd_community[, pd_tree$tip.label, drop = FALSE]
stopifnot(identical(colnames(pd_community), pd_tree$tip.label))
pd_taxa <- pd_tree$tip.label

# Precompute parent-edge lookup for efficient Faith's PD accumulation
pd_n_edges <- nrow(pd_tree$edge)
pd_edge_child <- pd_tree$edge[, 2]
pd_edge_parent <- pd_tree$edge[, 1]
pd_edge_lengths <- pd_tree$edge.length
pd_max_node <- max(pd_tree$edge)
pd_parent_edge <- integer(pd_max_node)
pd_parent_node <- integer(pd_max_node)
pd_parent_edge[pd_edge_child] <- seq_len(pd_n_edges)
pd_parent_node[pd_edge_child] <- pd_edge_parent
pd_tip_ids <- match(colnames(pd_community), pd_tree$tip.label)
if (anyNA(pd_tip_ids)) stop("Could not match community taxa to tree tips.")

# Named list required by the PD permutation test; function is defined
# here rather than relying on objects from a previous interactive session.
perm_sample_taxa <- lapply(seq_len(nrow(pd_community)), function(i) which(pd_community[i, ] > 0))
names(perm_sample_taxa) <- rownames(pd_community)
fast_pd_subset <- function(sample_ids) {
  taxa <- unique(unlist(perm_sample_taxa[sample_ids], use.names = FALSE))
  if (!length(taxa)) return(0)
  represented <- logical(pd_n_edges)
  total <- 0
  for (j in taxa) {
    node <- pd_tip_ids[j]
    repeat {
      if (node < 1L || node > pd_max_node) break
      e <- pd_parent_edge[node]
      if (e == 0L || represented[e]) break
      represented[e] <- TRUE
      total <- total + pd_edge_lengths[e]
      node <- pd_parent_node[node]
    }
  }
  total
}

# 7. FAITH'S PD ACCUMULATION CURVES ---------------------------
pd_curves <- run_curve_set("pd")

# 8. FAITH'S PD PAIRWISE PERMUTATION TESTS --------------------
pd_sample_metric <- function(ids) fast_pd_subset(ids)
pd_tests <- run_pairwise_tests("Faith's PD", pd_sample_metric, names(perm_sample_taxa))
write.csv(pd_tests, "Faith_PD_MaxCommonEffort_Pairwise_Tests.csv", row.names = FALSE)
openxlsx::write.xlsx(pd_tests, "Faith_PD_MaxCommonEffort_Pairwise_Tests.xlsx", overwrite = TRUE, rowNames = FALSE)

# 9. OUTPUT SUMMARY ------------------------------------------
cat("\nAll analyses completed. Files written to: ", getwd(), "\n", sep = "")
cat("\nNon-phylogenetic outputs:\n",
    "- Nonphylogenetic_Taxon_Accumulation_2024_2025_Both.pdf\n",
    "- Nonphylogenetic_Taxon_Accumulation_2024_2025_Both.png\n",
    "- Nonphylogenetic_Taxon_Accumulation_Summary_All_Years.csv\n",
    "- Nonphylogenetic_MaxCommonEffort_Pairwise_Tests.xlsx\n",
    "- Nonphylogenetic_MaxCommonEffort_Significant_Tests.xlsx\n", sep = "")
cat("\nFaith's PD outputs:\n",
    "- rZOTU_phylogenetic_tree.nwk\n",
    "- Faith_PD_Accumulation_2024_2025_Both.pdf\n",
    "- Faith_PD_Accumulation_2024_2025_Both.png\n",
    "- Faith_PD_Accumulation_Summary_All_Years.csv\n",
    "- Faith_PD_MaxCommonEffort_Pairwise_Tests.csv\n",
    "- Faith_PD_MaxCommonEffort_Pairwise_Tests.xlsx\n", sep = "")
