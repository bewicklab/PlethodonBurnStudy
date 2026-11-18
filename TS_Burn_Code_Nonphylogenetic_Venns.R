############################################################
# NON-PHYLOGENETIC GENUS CORE VENN ANALYSIS
#
# This follows the same workflow as the phylogenetic version:
#   - same metadata cleaning and sample filtering
#   - same control subtraction and rarefaction as original; then ASVs pooled to Genus
#   - 10 random subsampling trials
#   - 2025: 14 samples/site
#   - 2024: 10 samples/site
#   - combined: 14 samples/site in 2025 + 10/site in 2024
#   - average exclusive Venn-compartment values across trials
#
# Difference:
#   Core membership is based on genus prevalence (presence/absence)
#   within each site. A genus is core in a site when it occurs in
#   at least core_fraction of that site's samples.
#   Venn values are numbers of genera, not branch lengths.
############################################################

############################
# PACKAGES (retained from original workflow)
############################

library('sf')
library('ggplot2')
library('rnaturalearth')
library('rnaturalearthdata')
library(ggspatial)
library(rosm)
library(prettymapr)
library('maptiles')
library('ggmap')
library('tmap')
library('phytools')
library(qiime2R)
library(phyloseq)
library(dplyr)
library(readr)
library(lme4)
library("phyloseq")
library("ggplot2")
library("VennDiagram")
library(venneuler)
library(gplots)
library(vegan)
library("qiime2R")
library(devtools)
library(car)
library(tidyr)
library(dplyr)
library('stringr')
library('microshades')
library('speedyseq')
library('forcats')
library('cowplot')
library('knitr')
library('ggplot2')
library('PERMANOVA')
library('pairwiseAdonis')
library('ggsignif')
library('tidyverse')
library('ggpubr')
library('betapart')
library('stats')
library(ape)
library('geodist')
library(picante)
library(dplyr)
library(MASS)
library(lmerTest)
library("car")
library(cAIC4)
library("corrplot")
library(patchwork)
library("devtools")
library('holobiont')

suppressPackageStartupMessages({
  library(phyloseq)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(lefser)
  library(ggtree)
  library(ggtreeExtra)
  library(ggplot2)
  library(phangorn)
  library(ape)
})

############################################################
# NON-PHYLOGENETIC CORE GENUS VENN FUNCTION
############################################################

output_coreGenusVenn <- function(
    x,
    grouping,
    core_fraction = 0.5,
    ordered_groups = NULL,
    show_percentage = TRUE,
    decimal = 2,
    fill_color = c('red', 'orange', 'yellow', 'green',
                   'blue', 'purple', 'black'),
    fill_alpha = 0.5,
    stroke_color = 'black',
    stroke_alpha = 1,
    stroke_size = 1,
    stroke_linetype = "solid",
    set_name_color = "black",
    set_name_size = 6,
    text_color = "black",
    text_size = 4,
    return_values = FALSE
) {
  
  if (!inherits(x, "phyloseq")) {
    stop("'x' must be a phyloseq object.")
  }
  
  grouping <- as.character(grouping)
  
  if (length(grouping) != nsamples(x)) {
    stop("'grouping' must have one value for every sample in x.")
  }
  
  if (anyNA(grouping)) {
    stop("'grouping' contains NA values.")
  }
  
  if (length(core_fraction) != 1 ||
      is.na(core_fraction) ||
      core_fraction < 0 ||
      core_fraction > 1) {
    stop("'core_fraction' must be a single value between 0 and 1.")
  }
  
  group_id <- if (is.null(ordered_groups)) {
    unique(grouping)
  } else {
    as.character(ordered_groups)
  }
  
  if (length(group_id) < 2 || length(group_id) > 7) {
    stop("Venn diagrams require between 2 and 7 groups.")
  }
  
  if (!setequal(unique(grouping), group_id)) {
    stop("'ordered_groups' must contain exactly the groups in grouping.")
  }
  
  names(grouping) <- sample_names(x)
  
  otu_mat <- as(otu_table(x), "matrix")
  if (!taxa_are_rows(x)) {
    otu_mat <- t(otu_mat)
  }
  
  # Each row is now a genus (because the object was tax_glom'ed).
  # Presence means at least one read for that genus in that sample.
  presence_mat <- otu_mat > 0
  
  core_sets <- vector("list", length(group_id))
  names(core_sets) <- group_id
  
  for (i in seq_along(group_id)) {
    group_samples <- names(grouping)[grouping == group_id[i]]
    
    if (length(group_samples) == 0) {
      stop("No samples found for group: ", group_id[i])
    }
    
    group_presence <- presence_mat[
      ,
      group_samples,
      drop = FALSE
    ]
    
    prevalence <- rowMeans(group_presence)
    core_sets[[i]] <- rownames(presence_mat)[
      prevalence >= core_fraction
    ]
  }
  
  all_core_genera <- unique(unlist(core_sets, use.names = FALSE))
  
  if (length(all_core_genera) == 0) {
    stop("No core genera were identified in any group.")
  }
  
  membership <- as.data.frame(
    matrix(
      FALSE,
      nrow = length(all_core_genera),
      ncol = length(group_id),
      dimnames = list(all_core_genera, group_id)
    )
  )
  
  for (i in seq_along(group_id)) {
    membership[
      all_core_genera %in% core_sets[[i]],
      group_id[i]
    ] <- TRUE
  }
  
  membership <- membership[rowSums(membership) > 0, , drop = FALSE]
  
  combos <- list()
  compartment_counts <- numeric(0)
  
  for (i in seq_along(group_id)) {
    ff <- combn(seq_along(group_id), i)
    
    for (j in seq_len(ncol(ff))) {
      included <- ff[, j]
      outside <- setdiff(seq_along(group_id), included)
      
      in_all_included <- apply(
        membership[, included, drop = FALSE],
        1,
        all
      )
      
      in_any_outside <- if (length(outside) > 0) {
        apply(
          membership[, outside, drop = FALSE],
          1,
          any
        )
      } else {
        rep(FALSE, nrow(membership))
      }
      
      exclusive_genera <- rownames(membership)[
        in_all_included & !in_any_outside
      ]
      
      combos[[length(combos) + 1]] <- included
      compartment_counts <- c(
        compartment_counts,
        length(exclusive_genera)
      )
    }
  }
  
  compartment_values <- data.frame(
    Compartment = vapply(
      combos,
      function(idx) paste(group_id[idx], collapse = " & "),
      character(1)
    ),
    Genus_count = compartment_counts,
    stringsAsFactors = FALSE
  )
  
  total_core_genera <- sum(compartment_values$Genus_count)
  
  if (total_core_genera == 0) {
    stop("There are no core genera in any group.")
  }
  
  compartment_values$Percentage <-
    100 * compartment_values$Genus_count / total_core_genera
  
  compartment_values$Genus_count <-
    round(compartment_values$Genus_count, decimal)
  
  compartment_values$Percentage <-
    round(compartment_values$Percentage, decimal)
  
  venn_plot <- ggvenn2(
    membership,
    columns = group_id,
    show_elements = FALSE,
    show_percentage = show_percentage,
    digits = decimal,
    fill_color = fill_color,
    fill_alpha = fill_alpha,
    stroke_color = stroke_color,
    stroke_alpha = stroke_alpha,
    stroke_size = stroke_size,
    stroke_linetype = stroke_linetype,
    set_name_color = set_name_color,
    set_name_size = set_name_size,
    text_color = text_color,
    text_size = text_size
  ) + ggplot2::coord_cartesian(clip = "off")
  
  attr(venn_plot, "compartment_values") <- compartment_values
  attr(venn_plot, "core_sets") <- core_sets
  
  if (return_values) {
    return(list(
      plot = venn_plot,
      compartment_values = compartment_values,
      core_sets = core_sets,
      membership = membership
    ))
  }
  
  return(venn_plot)
}

############################################################
# HELPER: AVERAGE TRIAL VALUES AND PLOT AN AVERAGE VENN
############################################################

average_genus_trials <- function(trial_values,
                                 compartments,
                                 fill_colors,
                                 scale_factor = 100,
                                 set_name_size = 6,
                                 text_size = 4.8) {
  
  all_trial_values <- do.call(rbind, trial_values)
  
  # Ensure every trial contributes a value for every compartment.
  complete_values <- merge(
    expand.grid(
      Trial = sort(unique(all_trial_values$Trial)),
      Compartment = compartments,
      stringsAsFactors = FALSE
    ),
    all_trial_values,
    by = c("Trial", "Compartment"),
    all.x = TRUE,
    sort = FALSE
  )
  
  # A missing compartment in a trial means zero genera in that compartment.
  complete_values$Genus_count[is.na(complete_values$Genus_count)] <- 0
  
  average_values <- aggregate(
    Genus_count ~ Compartment,
    data = complete_values,
    FUN = mean
  )
  
  sd_values <- aggregate(
    Genus_count ~ Compartment,
    data = complete_values,
    FUN = sd
  )
  
  names(sd_values)[2] <- "SD"
  
  average_values <- merge(
    average_values,
    sd_values,
    by = "Compartment",
    sort = FALSE
  )
  
  average_values$Percentage <-
    100 * average_values$Genus_count /
    sum(average_values$Genus_count)
  
  average_values$Genus_count <- round(average_values$Genus_count, 6)
  average_values$SD <- round(average_values$SD, 6)
  average_values$Percentage <- round(average_values$Percentage, 2)
  
  average_values$Compartment <- factor(
    average_values$Compartment,
    levels = compartments
  )
  
  average_values <- average_values[
    order(average_values$Compartment),
  ]
  
  # Convert mean genus counts to integer membership rows for plotting.
  # Multiplying by 100 retains two decimal places in the averages.
  mean_counts <- setNames(
    average_values$Genus_count,
    as.character(average_values$Compartment)
  )
  mean_counts <- mean_counts[compartments]
  mean_counts[is.na(mean_counts)] <- 0
  
  plot_counts <- round(mean_counts * scale_factor)
  
  final_avg <- as.data.frame(
    matrix(
      FALSE,
      nrow = sum(plot_counts),
      ncol = length(c("Burn", "Control_2", "Control_1")),
      dimnames = list(
        NULL,
        c("Burn", "Control_2", "Control_1")
      )
    )
  )
  
  counter <- 1L
  
  for (i in seq_along(compartments)) {
    n <- plot_counts[i]
    
    if (n > 0) {
      idx <- seq.int(counter, length.out = n)
      
      groups <- strsplit(
        compartments[i],
        " & ",
        fixed = TRUE
      )[[1]]
      
      final_avg[idx, groups] <- TRUE
      counter <- counter + n
    }
  }
  
  if (nrow(final_avg) == 0) {
    stop("No core genera were available to plot.")
  }
  
  stopifnot(all(rowSums(final_avg) > 0))
  
  p_average <- ggvenn2(
    final_avg,
    columns = NULL,
    show_elements = FALSE,
    show_percentage = FALSE,
    digits = 2,
    fill_color = fill_colors,
    fill_alpha = 0.5,
    stroke_color = "black",
    stroke_alpha = 1,
    stroke_size = 1,
    stroke_linetype = "solid",
    set_name_color = "black",
    set_name_size = set_name_size,
    text_color = "black",
    text_size = text_size
  ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica")
    )
  
  list(
    all_trial_values = all_trial_values,
    average_values = average_values,
    plot = p_average
  )
}


############################################################
# LOAD & CLEAN METADATA
############################################################

setwd("ZTL_Sequence_Data_All/")

meta <- read.csv(
  "Coweeta_TEWL_Seasons_Sites.csv",
  stringsAsFactors = FALSE
)

meta <- meta %>%
  filter(QIIME_ID_aligned != "") %>%
  distinct(QIIME_ID_aligned, .keep_all = TRUE)

rownames(meta) <- meta$QIIME_ID_aligned
rownames(meta) <- tolower(rownames(meta))
rownames(meta) <- gsub("-", ".", rownames(meta))


############################################################
# IMPORT QIIME DATA
############################################################

# The tree is retained here to keep the import step consistent
# with the original script. The Genus Venn analysis below does not
# use the tree.
Z_physeq <- qza_to_phyloseq(
  features = "table.qza",
  taxonomy = "taxonomy.qza",
  tree = "tree.qza"
)


############################################################
# CLEAN SAMPLE NAMES AND MATCH METADATA
############################################################

sample_names(Z_physeq) <- tolower(sample_names(Z_physeq))
sample_names(Z_physeq) <- gsub("-", ".", sample_names(Z_physeq))

common <- intersect(sample_names(Z_physeq), rownames(meta))

Z_physeq <- prune_samples(common, Z_physeq)
meta <- meta[common, ]

sample_data(Z_physeq) <- sample_data(meta)

# Sanity checks
print(all(sample_names(Z_physeq) == rownames(sample_data(Z_physeq))))
print(setdiff(taxa_names(Z_physeq), phy_tree(Z_physeq)$tip.label))


############################################################
# REDUCE DATA SIZE
############################################################

Z_physeq <- prune_taxa(taxa_sums(Z_physeq) > 10, Z_physeq)


############################################################
# CONTROL SUBTRACTION (SAME AS ORIGINAL PHYLOGENETIC SCRIPT)
############################################################

control_labels <- c("Blank")

control_samples <- sample_names(Z_physeq)[
  sample_data(Z_physeq)$Site %in% control_labels
]

if (length(control_samples) == 0) {
  stop("No control samples detected.")
}

otu_mat <- as(otu_table(Z_physeq), "matrix")

if (!taxa_are_rows(Z_physeq)) {
  otu_mat <- t(otu_mat)
}

control_mat <- otu_mat[, control_samples, drop = FALSE]

max_control_reads <- apply(control_mat, 1, max, na.rm = TRUE)
max_control_reads[is.infinite(max_control_reads)] <- 0

corrected_mat <- sweep(otu_mat, 1, max_control_reads, FUN = "-")
corrected_mat[corrected_mat < 0] <- 0

otu_table(Z_physeq) <- otu_table(corrected_mat, taxa_are_rows = TRUE)


############################################################
# REMOVE CONTROL SAMPLES
############################################################

Z_physeq <- subset_samples(
  Z_physeq,
  !(Site %in% control_labels)
)

Z_physeq <- prune_taxa(taxa_sums(Z_physeq) > 0, Z_physeq)


############################################################
# REMOVE UNWANTED TAXA
############################################################

Z_physeq <- subset_taxa(
  Z_physeq,
  !(Kingdom %in% c("d__Eukaryota")) &
    !(Phylum %in% c(
      "Chytridiomycota", "Ascomycota", "Chlorophyta",
      "Apicomplexa", "Ciliophora", "Phragmoplastophyta",
      "Zoopagomycota"
    )) &
    Family != "Chloroplast" &
    Genus != "Mitochondria"
)

Z_physeq <- prune_taxa(taxa_sums(Z_physeq) > 0, Z_physeq)


############################################################
# REMOVE NON-BURN STUDY SAMPLES
############################################################

Z_physeq <- prune_samples(sample_data(Z_physeq)$Site != 'Control', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$Site != '2750', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$Site != '2850', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$Site != 'EF', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$Site != 'PN', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$Site != 'Positive_control', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$title != 'TS.31.5', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$title != 'TS.31.6', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$title != 'TS.31.7', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$title != 'TS.31.8', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$title != 'TS.31.9', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$title != 'TS.31.10', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$title != 'TS.31.11', Z_physeq)
Z_physeq <- prune_samples(
  sample_data(Z_physeq)$QIIME_ID_aligned != 'C1-25-27-NoBurn',
  Z_physeq
)
Z_physeq <- prune_samples(
  sample_data(Z_physeq)$Season != 'Summer_2024',
  Z_physeq
)
Z_physeq <- prune_samples(
  sample_data(Z_physeq)$Season != 'Summer_2025',
  Z_physeq
)
Z_physeq <- prune_samples(
  sample_data(Z_physeq)$Season != 'Fall_2024',
  Z_physeq
)
Z_physeq <- prune_samples(
  !is.na(sample_data(Z_physeq)$lat_7_Collection_locatio),
  Z_physeq
)


############################################################
# RAREFY (SAME AS ORIGINAL PHYLOGENETIC SCRIPT)
############################################################

set.seed(123)

rZOTU <- rarefy_even_depth(
  Z_physeq,
  sample.size = 12000,
  rngseed = 123,
  replace = FALSE
)


############################################################
# POOL RAREFIED ASVs TO GENUS
############################################################

# Exclude taxa without an informative genus assignment rather than
# combining all unclassified ASVs into one artificial "unclassified" genus.
genus_labels <- as.character(tax_table(rZOTU)[, "Genus"])
genus_labels_trimmed <- trimws(genus_labels)

valid_genus <- !is.na(genus_labels_trimmed) &
  genus_labels_trimmed != "" &
  !tolower(genus_labels_trimmed) %in% c(
    "na", "unassigned", "uncultured", "unknown",
    "unclassified", "g__", "genus"
  )

rZOTU <- prune_taxa(valid_genus, rZOTU)

if (ntaxa(rZOTU) == 0) {
  stop("No taxa with usable Genus assignments remain after filtering.")
}

# Pool ASVs that share the same genus-level taxonomy.
rZOTU <- tax_glom(
  rZOTU,
  taxrank = "Genus",
  NArm = TRUE
)

rZOTU <- prune_taxa(taxa_sums(rZOTU) > 0, rZOTU)

# Give taxa genus-based names for readable membership outputs.
genus_names <- as.character(tax_table(rZOTU)[, "Genus"])
taxa_names(rZOTU) <- make.unique(genus_names)

message("Number of genus-level taxa after pooling: ", ntaxa(rZOTU))


# Rename sites
sample_data(rZOTU)$Site[
  sample_data(rZOTU)$Site == '31_burn'
] <- 'Burn'

sample_data(rZOTU)$Site[
  sample_data(rZOTU)$Site == '32_no_burn'
] <- 'Control_2'

sample_data(rZOTU)$Site[
  sample_data(rZOTU)$Site == '2950'
] <- 'Control_1'

SitexYear <- c()
Year <- c()

for (k in seq_along(sample_data(rZOTU)$Site)) {
  SitexYear <- c(
    SitexYear,
    paste0(
      sample_data(rZOTU)$Site[k], '_',
      strsplit(sample_data(rZOTU)$Season[k], '_')[[1]][2]
    )
  )
  
  Year <- c(
    Year,
    strsplit(sample_data(rZOTU)$Season[k], '_')[[1]][2]
  )
}

sample_data(rZOTU)$SitexYear <- SitexYear
sample_data(rZOTU)$Year <- Year

sample_data(rZOTU)$BurnStatus <- ifelse(
  sample_data(rZOTU)$Site == "Burn",
  "Burn",
  "Control"
)

sample_data(rZOTU)$BurnStatus <- factor(
  sample_data(rZOTU)$BurnStatus,
  levels = c("Control", "Burn")
)


############################################################
# SHARED SETTINGS
############################################################

n_trials <- 10
core_fraction <- 0.9

set.seed(12345)

ordered_groups <- c("Burn", "Control_2", "Control_1")
fill_colors <- c("red", "green", "forestgreen")

compartments <- c(
  "Burn",
  "Control_2",
  "Control_1",
  "Burn & Control_2",
  "Burn & Control_1",
  "Control_2 & Control_1",
  "Burn & Control_2 & Control_1"
)


############################################################
# 2025 ANALYSIS: 14 SAMPLES PER SITE
############################################################

n_per_site <- 14

rZOTU_2025 <- prune_samples(
  sample_data(rZOTU)$Season != "Spring_2024",
  rZOTU
)

samples_2025_by_site <- split(
  sample_names(rZOTU_2025),
  sample_data(rZOTU_2025)$Site
)

print(lengths(samples_2025_by_site))

if (any(lengths(samples_2025_by_site[ordered_groups]) < n_per_site)) {
  stop("At least one 2025 site has fewer than 14 available samples.")
}

trial_values_2025 <- vector("list", n_trials)

for (trial in seq_len(n_trials)) {
  
  message("2025: running trial ", trial, " of ", n_trials)
  
  keep <- unlist(
    lapply(
      samples_2025_by_site[ordered_groups],
      function(samples) sample(samples, n_per_site)
    ),
    use.names = FALSE
  )
  
  rZOTU_trial <- prune_samples(keep, rZOTU_2025)
  
  trial_plot <- output_coreGenusVenn(
    rZOTU_trial,
    sample_data(rZOTU_trial)$Site,
    core_fraction = core_fraction,
    ordered_groups = ordered_groups,
    fill_color = fill_colors,
    show_percentage = FALSE
  )
  
  values <- attr(trial_plot, "compartment_values")
  
  if (is.null(values)) {
    stop("2025 trial ", trial, " did not return compartment_values.")
  }
  
  values$Trial <- trial
  trial_values_2025[[trial]] <- values
}

result_2025 <- average_genus_trials(
  trial_values_2025,
  compartments,
  fill_colors
)

all_trial_values_2025 <- result_2025$all_trial_values
average_values_2025 <- result_2025$average_values
p_average_2025 <- result_2025$plot

print(average_values_2025)

write.csv(
  all_trial_values_2025,
  "2025_nonphylo_genus_venn_10_trials.csv",
  row.names = FALSE
)

write.csv(
  average_values_2025,
  "2025_nonphylo_genus_venn_averages.csv",
  row.names = FALSE
)


############################################################
# 2024 ANALYSIS: 10 SAMPLES PER SITE
############################################################

n_per_site <- 10

# Preserve the original year-selection logic:
# 2024 includes samples other than Spring_2025.
rZOTU_2024 <- prune_samples(
  sample_data(rZOTU)$Season != "Spring_2025",
  rZOTU
)

samples_2024_by_site <- split(
  sample_names(rZOTU_2024),
  sample_data(rZOTU_2024)$Site
)

print(lengths(samples_2024_by_site))

if (any(lengths(samples_2024_by_site[ordered_groups]) < n_per_site)) {
  stop("At least one 2024 site has fewer than 10 available samples.")
}

trial_values_2024 <- vector("list", n_trials)

for (trial in seq_len(n_trials)) {
  
  message("2024: running trial ", trial, " of ", n_trials)
  
  keep <- unlist(
    lapply(
      samples_2024_by_site[ordered_groups],
      function(samples) sample(samples, n_per_site)
    ),
    use.names = FALSE
  )
  
  rZOTU_trial <- prune_samples(keep, rZOTU_2024)
  
  trial_plot <- output_coreGenusVenn(
    rZOTU_trial,
    sample_data(rZOTU_trial)$Site,
    core_fraction = core_fraction,
    ordered_groups = ordered_groups,
    fill_color = fill_colors,
    show_percentage = FALSE
  )
  
  values <- attr(trial_plot, "compartment_values")
  
  if (is.null(values)) {
    stop("2024 trial ", trial, " did not return compartment_values.")
  }
  
  values$Trial <- trial
  trial_values_2024[[trial]] <- values
}

result_2024 <- average_genus_trials(
  trial_values_2024,
  compartments,
  fill_colors
)

all_trial_values_2024 <- result_2024$all_trial_values
average_values_2024 <- result_2024$average_values
p_average_2024 <- result_2024$plot

print(average_values_2024)

write.csv(
  all_trial_values_2024,
  "2024_nonphylo_genus_venn_10_trials.csv",
  row.names = FALSE
)

write.csv(
  average_values_2024,
  "2024_nonphylo_genus_venn_averages.csv",
  row.names = FALSE
)


############################################################
# COMBINED 2024 + 2025 ANALYSIS
# 14 samples/site from 2025 + 10 samples/site from 2024
############################################################

samples_2025_by_site <- split(
  sample_names(rZOTU_2025),
  sample_data(rZOTU_2025)$Site
)

samples_2024_by_site <- split(
  sample_names(rZOTU_2024),
  sample_data(rZOTU_2024)$Site
)

print(lengths(samples_2024_by_site))
print(lengths(samples_2025_by_site))

if (any(lengths(samples_2024_by_site[ordered_groups]) < 10)) {
  stop("At least one 2024 site has fewer than 10 samples.")
}

if (any(lengths(samples_2025_by_site[ordered_groups]) < 14)) {
  stop("At least one 2025 site has fewer than 14 samples.")
}

trial_values_2024_2025 <- vector("list", n_trials)

for (trial in seq_len(n_trials)) {
  
  message("Combined 2024-2025: running trial ", trial, " of ", n_trials)
  
  keep_2025 <- unlist(
    lapply(
      samples_2025_by_site[ordered_groups],
      function(samples) sample(samples, 14)
    ),
    use.names = FALSE
  )
  
  keep_2024 <- unlist(
    lapply(
      samples_2024_by_site[ordered_groups],
      function(samples) sample(samples, 10)
    ),
    use.names = FALSE
  )
  
  keep <- c(keep_2025, keep_2024)
  
  rZOTU_trial <- prune_samples(keep, rZOTU)
  
  trial_plot <- output_coreGenusVenn(
    rZOTU_trial,
    sample_data(rZOTU_trial)$Site,
    core_fraction = core_fraction,
    ordered_groups = ordered_groups,
    fill_color = fill_colors,
    show_percentage = FALSE
  )
  
  values <- attr(trial_plot, "compartment_values")
  
  if (is.null(values)) {
    stop(
      "Combined trial ", trial,
      " did not return compartment_values."
    )
  }
  
  values$Trial <- trial
  trial_values_2024_2025[[trial]] <- values
}

result_2024_2025 <- average_genus_trials(
  trial_values_2024_2025,
  compartments,
  fill_colors,
  set_name_size = 6,
  text_size = 4.8
)

all_trial_values_2024_2025 <- result_2024_2025$all_trial_values
average_values_2024_2025 <- result_2024_2025$average_values
p_average_2024_2025 <- result_2024_2025$plot

print(average_values_2024_2025)

write.csv(
  all_trial_values_2024_2025,
  "2024_2025_nonphylo_genus_venn_10_trials.csv",
  row.names = FALSE
)

write.csv(
  average_values_2024_2025,
  "2024_2025_nonphylo_genus_venn_averages.csv",
  row.names = FALSE
)


############################################################
# THREE-PANEL FIGURE AND PDF EXPORT
# Layout:
#   2024       top left
#   2025       top right
#   2024+2025  bottom left
#   bottom right intentionally blank
############################################################

# Apply identical titles and typography to each panel.
panel_theme <- ggplot2::theme(
  text = ggplot2::element_text(family = "Helvetica"),
  plot.title = ggplot2::element_text(
    family = "Helvetica",
    size = 14,
    face = "bold",
    hjust = 0.5,
    margin = ggplot2::margin(b = 6)
  ),
  plot.margin = ggplot2::margin(10, 10, 10, 10)
)

p_2024_panel <- p_average_2024 +
  ggplot2::ggtitle("2024") +
  panel_theme

p_2025_panel <- p_average_2025 +
  ggplot2::ggtitle("2025") +
  panel_theme

p_2024_2025_panel <- p_average_2024_2025 +
  ggplot2::ggtitle("2024 + 2025") +
  panel_theme

three_panel_figure <- cowplot::plot_grid(
  p_2024_panel,
  p_2025_panel,
  p_2024_2025_panel,
  NULL,
  ncol = 2,
  nrow = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

# Save a landscape, letter-size PDF in the current working directory.
ggplot2::ggsave(
  filename = "Genus_Venn_Three_Panel.pdf",
  plot = three_panel_figure,
  device = function(filename, width, height, ...) {
    grDevices::pdf(
      file = filename,
      width = width,
      height = height,
      family = "Helvetica",
      useDingbats = FALSE,
      ...
    )
  },
  width = 11,
  height = 8.5,
  units = "in"
)

print(three_panel_figure)

message(
  "Three-panel PDF saved to: ",
  file.path(getwd(), "Genus_Venn_Three_Panel.pdf")
)

############################################################
# END
#
# Plot order:
#   1. p_average_2025
#   2. p_average_2024
#   3. p_average_2024_2025
#
# Tables:
#   average_values_2025
#   average_values_2024
#   average_values_2024_2025
############################################################
