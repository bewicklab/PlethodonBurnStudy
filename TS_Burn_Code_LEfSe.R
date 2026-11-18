

library(qiime2R)
library(phyloseq)
library(dplyr)
library(readr)
library(lme4)
library("phyloseq")
library("ggplot2")  #Used for Graphics
library("VennDiagram") #Used for Venn Diagram
library(venneuler) #Used for Venn Diagram
library(gplots) #Used for Graphics
#library(limma)
library(vegan)
#library("metagMisc")
library("qiime2R")
library(devtools)
library(car)
#library(metagenomeSeq)
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

########################
# LOAD & CLEAN METADATA
########################
setwd("ZTL_Sequence_Data_All/")
meta <- read.csv("Coweeta_TEWL_Seasons_Sites.csv", stringsAsFactors = FALSE)

# Remove empty IDs
meta <- meta %>%
  filter(QIIME_ID_aligned != "")

# Remove duplicate IDs
meta <- meta %>%
  distinct(QIIME_ID_aligned, .keep_all = TRUE)

# Set rownames
rownames(meta) <- meta$QIIME_ID_aligned

# Standardize names (VERY IMPORTANT)
rownames(meta) <- tolower(rownames(meta))
rownames(meta) <- gsub("-", ".", rownames(meta))

########################
# IMPORT QIIME DATA (NO METADATA YET)
########################

Z_physeq <- qza_to_phyloseq(
  features = "table.qza",
  taxonomy = "taxonomy.qza",
  tree = "tree.qza"
)



########################
# CLEAN SAMPLE NAMES IN PHYLOSEQ
########################

sample_names(Z_physeq) <- tolower(sample_names(Z_physeq))
sample_names(Z_physeq) <- gsub("-", ".", sample_names(Z_physeq))

########################
# MATCH SAMPLES BETWEEN DATA & METADATA
########################

common <- intersect(sample_names(Z_physeq), rownames(meta))

Z_physeq <- prune_samples(common, Z_physeq)
meta <- meta[common, ]

# Attach metadata
sample_data(Z_physeq) <- sample_data(meta)

########################
# SANITY CHECKS (DO NOT SKIP)
########################

# Should be TRUE
print(all(sample_names(Z_physeq) == rownames(sample_data(Z_physeq))))

# Should be empty
print(setdiff(taxa_names(Z_physeq), phy_tree(Z_physeq)$tip.label))

########################
# OPTIONAL: REDUCE DATA SIZE (HIGHLY RECOMMENDED)
########################

Z_physeq <- prune_taxa(taxa_sums(Z_physeq) > 10, Z_physeq)

########################
# CONTROL SUBTRACTION
########################

control_labels <- c("Blank")

control_samples <- sample_names(Z_physeq)[
  sample_data(Z_physeq)$Site %in% control_labels
]

if(length(control_samples) == 0){
  stop("No control samples detected.")
}

otu_mat <- as(otu_table(Z_physeq), "matrix")

if(!taxa_are_rows(Z_physeq)){
  otu_mat <- t(otu_mat)
}

control_mat <- otu_mat[, control_samples, drop = FALSE]

max_control_reads <- apply(control_mat, 1, max, na.rm = TRUE)
max_control_reads[is.infinite(max_control_reads)] <- 0

corrected_mat <- sweep(otu_mat, 1, max_control_reads, FUN = "-")
corrected_mat[corrected_mat < 0] <- 0

otu_table(Z_physeq) <- otu_table(corrected_mat, taxa_are_rows = TRUE)

########################
# REMOVE CONTROL SAMPLES
########################

Z_physeq <- subset_samples(
  Z_physeq,
  !(Site %in% control_labels)
)

Z_physeq <- prune_taxa(taxa_sums(Z_physeq) > 0, Z_physeq)

########################
# REMOVE UNWANTED TAXA
########################

Z_physeq <- subset_taxa(
  Z_physeq,
  !(Kingdom %in% c(
    "d__Eukaryota"
  )) &
    !(Phylum %in% c(
      "Chytridiomycota","Ascomycota","Chlorophyta",
      "Apicomplexa","Ciliophora","Phragmoplastophyta","Zoopagomycota"
    )) &
    Family != "Chloroplast" &
    Genus != "Mitochondria"
)

Z_physeq <- prune_taxa(taxa_sums(Z_physeq) > 0, Z_physeq)

########################
# FILTER LOW-DEPTH SAMPLES
########################

#Z_physeq <- prune_samples(sample_sums(Z_physeq) >= 2000, Z_physeq)


########################
# REMOVE NON-BURN STUDY SAMPLES
########################

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
Z_physeq <- prune_samples(sample_data(Z_physeq)$Season != 'Summer_2024', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$Season != 'Summer_2025', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$Season != 'Fall_2024', Z_physeq)
Z_physeq <- prune_samples(!(is.na(sample_data(Z_physeq)$lat_7_Collection_locatio)), Z_physeq)

########################
# RAREFY
########################

set.seed(123)

rZOTU <- rarefy_even_depth(
  Z_physeq,
  sample.size = 12000,
  rngseed = 123,
  replace = FALSE
)

########################
# CORRECT TAXON NAMES
########################


f1<-which(tax_table(rZOTU)[, "Genus"] == 'uncultured')
f2<-which(tax_table(rZOTU)[, "Family"] == 'uncultured')
f3<-which(tax_table(rZOTU)[, "Order"] == 'uncultured')
f4<-which(tax_table(rZOTU)[, "Class"] == 'uncultured')
f5<-which(tax_table(rZOTU)[, "Phylum"] == 'uncultured')

f1i<-which(!(f1 %in% f2))

suppressMessages({
  for (k in 1:length(f1i)) {
    tax_table(rZOTU)[f1[f1i[k]], 'Genus'] <-
      paste0('uncultured_', tax_table(rZOTU)[f1[f1i[k]], 'Family'])
  }
})


f2i<-which(!(f2 %in% f3))

suppressMessages(suppressWarnings(
  for (k in 1:length(f2i)) {
    tax_table(rZOTU)[f2[f2i[k]], 'Genus'] <-
      paste0('uncultured_', tax_table(rZOTU)[f2[f2i[k]], 'Order'])
  }
))

f3i<-which(!(f3 %in% f4))

suppressMessages({
  for (k in seq_along(f3i)) {
    tax_table(rZOTU)[f3[f3i[k]], 'Genus'] <-
      paste0('uncultured_', tax_table(rZOTU)[f3[f3i[k]], 'Class'])
  }
})

f4i<-which(!(f4 %in% f5))

suppressMessages({
  for (k in seq_along(f4i)) {
    tax_table(rZOTU)[f4[f4i[k]], 'Genus'] <-
      paste0('uncultured_', tax_table(rZOTU)[f4[f4i[k]], 'Phylum'])
  }
})

########################
# CORRECT SITE NAMES
########################

sample_data(rZOTU)$Site[which(sample_data(rZOTU)$Site == '31_burn')]<-'Burn'
sample_data(rZOTU)$Site[which(sample_data(rZOTU)$Site == '32_no_burn')]<-'Control_2'
sample_data(rZOTU)$Site[which(sample_data(rZOTU)$Site == '2950')]<-'Control_1'

########################
# ADD A SITExYEAR COLUMN
########################

SitexYear<-c()
for (k in 1:length(sample_data(rZOTU)$Site)){
  SitexYear<-c(SitexYear,paste0(sample_data(rZOTU)$Site[k],'_',strsplit(sample_data(rZOTU)$Season[k],'_')[[1]][2])) 
}

sample_data(rZOTU)$SitexYear<-SitexYear

table(sample_data(rZOTU)$BurnStatus)


rGenus <- tax_glom(
  rZOTU,
  taxrank = "Genus",
  NArm = FALSE
)


rGenus <- filter_taxa(
  rGenus,
  function(x) sum(x > 0) >= 3,
  prune = TRUE
)


ntaxa(rGenus)

library(phyloseq)
library(SummarizedExperiment)
library(S4Vectors)
library(lefser)

############################
# ABUNDANCE MATRIX
############################

otu <- as(otu_table(rGenus), "matrix")

if (!taxa_are_rows(rGenus)) {
  otu <- t(otu)
}

############################
# METADATA
############################

meta <- as.data.frame(sample_data(rGenus))

############################
# TAXONOMY
############################

tax <- as.data.frame(tax_table(rGenus))

############################
# CREATE LEfSe DATA
# samples = rows
# genera = columns
############################

lefse_data <- as.data.frame(t(otu))

############################
# BURN vs CONTROL
############################

site <- as.character(meta$Site)

lefse_data$BurnStatus <- ifelse(
  site == "Burn",
  "Burn",
  "Control"
)

lefse_data$BurnStatus <- factor(
  lefse_data$BurnStatus,
  levels = c("Control", "Burn")
)

############################
# CHECK GROUPS
############################

print(table(lefse_data$BurnStatus))

############################
# 5% PREVALENCE FILTER
############################

prev_threshold <- 0.05 * nrow(lefse_data)

taxa_cols <- setdiff(
  colnames(lefse_data),
  "BurnStatus"
)

keep <- colSums(
  lefse_data[, taxa_cols, drop = FALSE] > 0
) >= prev_threshold

lefse_data <- lefse_data[
  ,
  c(which(keep),
    which(colnames(lefse_data) == "BurnStatus")),
  drop = FALSE
]

############################
# CHECK FILTERING
############################

print(dim(lefse_data))

############################
# CREATE FEATURE MATRIX
# rows = genera
# columns = samples
############################

taxa_cols <- setdiff(
  colnames(lefse_data),
  "BurnStatus"
)

otu_lefse <- t(
  as.matrix(
    lefse_data[, taxa_cols, drop = FALSE]
  )
)

############################
# MAKE GENUS NAMES UNIQUE
#
# Prefix with "Genus__" so that
# lefser does not interpret one
# genus name as being nested
# inside another genus name.
############################

rownames(otu_lefse) <- paste0(
  "Genus__",
  rownames(otu_lefse)
)

############################
# CREATE SUMMARIZED EXPERIMENT
############################

se <- SummarizedExperiment(
  assays = list(
    counts = otu_lefse
  ),
  colData = DataFrame(
    BurnStatus = lefse_data$BurnStatus
  )
)

############################
# CHECK
############################

print(dim(se))

print(table(se$BurnStatus))

print(head(rownames(se)))

############################
# CONVERT TO RELATIVE
# ABUNDANCE
############################

se_rel <- relativeAb(se)

############################
# CHECK FOR DUPLICATE GENERA
############################

sum(duplicated(rownames(se_rel)))

############################
# REMOVE ZERO-VARIANCE GENERA
############################

x <- assay(se_rel)

keep_var <- apply(
  x,
  1,
  function(z) var(z) > 0
)

se_rel <- se_rel[keep_var, ]

############################
# CHECK FOR PERFECTLY
# COLLINEAR GENERA
############################

x <- assay(se_rel)

qr_x <- qr(t(x))

cat("Genera:", nrow(x), "\n")
cat("Samples:", ncol(x), "\n")
cat("Matrix rank:", qr_x$rank, "\n")

############################
# LEfSe WITH FDR CORRECTION
############################

lefse_result <- lefser(
  se_rel,
  classCol = "BurnStatus",
  lda.threshold = 2,
  method = "BH"
)

############################
# RESULTS
############################

lefse_result

############################
# MAP LEfSe FEATURES
# BACK TO GENUS NAMES
############################

lefse_df <- as.data.frame(lefse_result)

lefse_df$FeatureID <- sub(
  "^Genus__",
  "",
  lefse_df$features
)

tax_df <- as.data.frame(tax_table(rGenus))

tax_df$FeatureID <- rownames(tax_df)

lefse_tax <- merge(
  lefse_df,
  tax_df,
  by = "FeatureID",
  all.x = TRUE,
  sort = FALSE
)

############################
# ORDER BY LDA SCORE
############################

lefse_tax <- lefse_tax[
  order(-abs(lefse_tax$scores)),
]

############################
# LOOK AT RESULTS
############################

lefse_tax[
  ,
  c(
    "scores",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus"
  )
]


############################
# LEfSe RESULTS + ABUNDANCE
############################

# LEfSe results
lefse_df <- as.data.frame(lefse_result)

# Remove "Genus__" prefix to get original feature IDs
lefse_df$FeatureID <- sub(
  "^Genus__",
  "",
  lefse_df$features
)

# Relative abundance matrix
x <- assay(se_rel)

# Sample groups
group <- se_rel$BurnStatus

############################
# CALCULATE GROUP ABUNDANCES
############################

results <- lapply(seq_len(nrow(lefse_df)), function(i) {
  
  feature <- lefse_df$features[i]
  
  values <- x[feature, ]
  
  burn <- values[group == "Burn"]
  control <- values[group == "Control"]
  
  data.frame(
    Feature = feature,
    LDA = lefse_df$scores[i],
    
    Burn_Mean = mean(burn),
    Control_Mean = mean(control),
    
    Burn_Median = median(burn),
    Control_Median = median(control),
    
    Direction = ifelse(
      mean(burn) > mean(control),
      "Burn",
      "Control"
    ),
    
    stringsAsFactors = FALSE
  )
})

results <- do.call(rbind, results)

############################
# ADD TAXONOMY
############################

tax_df <- as.data.frame(tax_table(rGenus))

tax_df$FeatureID <- rownames(tax_df)

results$FeatureID <- sub(
  "^Genus__",
  "",
  results$Feature
)

results <- merge(
  results,
  tax_df,
  by = "FeatureID",
  all.x = TRUE,
  sort = FALSE
)

############################
# ORDER BY LDA EFFECT SIZE
############################

results <- results[
  order(-abs(results$LDA)),
]

############################
# SELECT COLUMNS
############################

results <- results[
  ,
  c(
    "LDA",
    "Direction",
    "Burn_Mean",
    "Control_Mean",
    "Burn_Median",
    "Control_Median",
    "Genus",
    "Family",
    "Phylum"
    
    
  )
]

############################
# VIEW RESULTS
############################

results


results_filtered <- results[
  (results$Burn_Median != 0 | results$Control_Median != 0) &
    abs(results$Burn_Median - results$Control_Median) /
    pmax(results$Burn_Median, results$Control_Median) > 0.20,
  ,
  drop = FALSE
]


results_filtered
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################









library(phyloseq)
library(ggtree)
library(ggplot2)
library(phangorn)

############################
# PHYLOGENETIC TREE
############################

tree <- phy_tree(rGenus)

############################
# TAXONOMY
############################

tax <- as.data.frame(tax_table(rGenus))

############################
# LEfSe TAXA
############################

burn_genera <- results_filtered$Genus[
  results_filtered$Direction == "Burn"
]

control_genera <- results_filtered$Genus[
  results_filtered$Direction == "Control"
]

############################
# TAXONOMY FOR TREE TIPS
############################

tree_tax <- data.frame(
  OTU = rownames(tax),
  Genus = as.character(tax$Genus),
  stringsAsFactors = FALSE
)

############################
# MATCH TAXONOMY TO TREE
############################

tree_tax <- tree_tax[
  match(tree$tip.label, tree_tax$OTU),
]

############################
# TIP GROUPS
############################

tree_tax$Group <- "Other"

tree_tax$Group[
  tree_tax$Genus %in% burn_genera
] <- "Burn"

tree_tax$Group[
  tree_tax$Genus %in% control_genera
] <- "Control"

############################
# TIP GROUP LOOKUP
############################

tip_group <- setNames(
  tree_tax$Group,
  tree_tax$OTU
)

############################
# CALCULATE GROUP FOR
# EACH TREE NODE
############################

node_group <- data.frame(
  node = seq_len(
    length(tree$tip.label) + tree$Nnode
  ),
  Group = "Other",
  stringsAsFactors = FALSE
)

for (node in seq_len(
  length(tree$tip.label) + tree$Nnode
)) {
  
  # Tip
  if (node <= length(tree$tip.label)) {
    
    group <- tip_group[
      tree$tip.label[node]
    ]
    
    if (!is.na(group)) {
      node_group$Group[
        node_group$node == node
      ] <- group
    }
    
  } else {
    
    # Internal node
    descendants <- phangorn::Descendants(
      tree,
      node,
      type = "tips"
    )[[1]]
    
    descendant_names <- tree$tip.label[
      descendants
    ]
    
    groups <- tip_group[
      descendant_names
    ]
    
    groups <- groups[
      !is.na(groups)
    ]
    
    if (
      length(groups) > 0 &&
      all(groups == "Burn")
    ) {
      
      node_group$Group[
        node_group$node == node
      ] <- "Burn"
      
    } else if (
      length(groups) > 0 &&
      all(groups == "Control")
    ) {
      
      node_group$Group[
        node_group$node == node
      ] <- "Control"
    }
  }
}

############################
# MAKE GGTREE
############################

p <- ggtree(tree)

############################
# ADD NODE GROUPS TO
# GGTREE DATA
############################

p$data$Group <- node_group$Group[
  match(
    p$data$node,
    node_group$node
  )
]

############################
# FAN TREE
############################

p <- ggtree(
  tree,
  layout = "fan",
  open.angle = 5,
  branch.length = "none",
  color = "grey70",
  linewidth = 0.25
)

############################
# ADD GROUP AND GENUS
############################

p$data$Group <- node_group$Group[
  match(p$data$node, node_group$node)
]

p$data$Genus <- tree_tax$Genus[
  match(p$data$label, tree_tax$OTU)
]

############################
# SIGNIFICANT GENERA
############################

sig_genera <- unique(results_filtered$Genus)

############################
# PLOT
############################

p <- p +
  
  
  # colored balls only on significant tips
  geom_tippoint(
    data = subset(p$data, isTip & Group != "Other"),
    aes(color = Group),
    size = 2
  ) +
  
  
  scale_color_manual(
    values = c(
      "Other" = "grey70",
      "Burn" = "red",
      "Control" = "blue"
    )
  ) +
  
  # scale_linewidth_manual(
  #   values = c(
  #     "Other" = 0.25,
  #     "Burn" = 0.7,
  #     "Control" = 0.7
  #   )
  # ) +
  
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  )

p


############################
# TAXONOMY DIRECTLY FROM PHYLOSEQ
############################

tax2 <- as.data.frame(tax_table(rGenus), stringsAsFactors = FALSE)
tree_tax2 <- data.frame(
  OTU = tree$tip.label,
  Phylum = as.character(
    tax2[match(tree$tip.label, rownames(tax2)), "Phylum"]
  ),
  stringsAsFactors = FALSE
)

############################
# CHECK
############################

table(is.na(tree_tax2$Phylum))

head(tree_tax2)

#############################################################
# CORE PHYLOGENETIC CLADES
#
# A core clade must be at least 90% one phylum.
#
# We choose the largest such clade for each phylum.
#############################################################

core_phyla <- c(
  "Actinobacteriota",
  "Acidobacteriota",
  "Chloroflexi",
  "Bacteroidota",
  "Planctomycetota",
  "Cyanobacteria",
  "Gemmatimonadota",
  "Armatimonadota",
  "Proteobacteria",
  "Firmicutes",
  "Verrucomicrobiota",
  "Bdellovibrionota+Myxococcota"
)


#############################################################
# TAXONOMY
#############################################################

tax_all <- as.data.frame(
  tax_table(rGenus),
  stringsAsFactors = FALSE
)

tax_all$FeatureID <- rownames(tax_all)

# Pool these two phyla for the core-phyla tree band only.
tax_all$Phylum <- as.character(tax_all$Phylum)
tax_all$Phylum[
  tax_all$Phylum %in% c("Bdellovibrionota", "Myxococcota")
] <- "Bdellovibrionota+Myxococcota"


tree_tax <- data.frame(
  OTU = tree$tip.label,
  
  Phylum = as.character(
    tax_all$Phylum[
      match(
        tree$tip.label,
        tax_all$FeatureID
      )
    ]
  ),
  
  stringsAsFactors = FALSE
)


#############################################################
# FIND CANDIDATE NODES
#############################################################

n_tips <- length(tree$tip.label)

n_nodes <- n_tips + tree$Nnode

candidate_clades <- data.frame(
  node = integer(0),
  Phylum = character(0),
  n_tips = integer(0),
  n_target = integer(0),
  n_other = integer(0),
  proportion_target = numeric(0),
  stringsAsFactors = FALSE
)


for (node in (n_tips + 1):n_nodes) {
  
  tips <- phangorn::Descendants(
    tree,
    node,
    type = "tips"
  )[[1]]
  
  if (length(tips) < 2) {
    next
  }
  
  phyla <- tree_tax$Phylum[tips]
  
  for (ph in core_phyla) {
    
    n_target <- sum(
      phyla == ph,
      na.rm = TRUE
    )
    
    n_other <- sum(
      phyla != ph,
      na.rm = TRUE
    )
    
    n_missing <- sum(
      is.na(phyla)
    )
    
    n_other_total <- n_other + n_missing
    
    proportion_target <-
      n_target / length(tips)
    
    #########################################################
    # CORE = >=90% TARGET PHYLUM
    #########################################################
    
    if (
      n_target >= 5 &&
      proportion_target >= 0.70
    ) {
      
      candidate_clades <- rbind(
        candidate_clades,
        data.frame(
          node = node,
          Phylum = ph,
          n_tips = length(tips),
          n_target = n_target,
          n_other = n_other_total,
          proportion_target =
            proportion_target,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}


#############################################################
# SELECT LARGEST TARGET CLade FOR EACH PHYLUM
#############################################################

if (nrow(candidate_clades) > 0) {
  core_clades <- candidate_clades[
    order(
      candidate_clades$Phylum,
      -candidate_clades$n_target,
      -candidate_clades$proportion_target
    ),
    ,
    drop = FALSE
  ]
  
  core_clades <- core_clades[
    !duplicated(core_clades$Phylum),
    ,
    drop = FALSE
  ]
  
  core_clades <- core_clades[
    order(match(core_clades$Phylum, core_phyla)),
    ,
    drop = FALSE
  ]
} else {
  core_clades <- candidate_clades
  message("No core clades met the selection criteria.")
}


#############################################################
# SHOW WHAT WAS SELECTED
#############################################################

print(core_clades)

#############################################################
# BUILD CORE PHYLA BAND
#############################################################

core_band <- data.frame(
  label = tree$tip.label,
  Phylum = NA_character_,
  stringsAsFactors = FALSE
)


if (nrow(core_clades) > 0) {
  for (i in seq_len(nrow(core_clades))) {
    tips <- phangorn::Descendants(
      tree,
      core_clades$node[i],
      type = "tips"
    )[[1]]
    
    core_band$Phylum[tips] <- core_clades$Phylum[i]
  }
}


#############################################################
# CHECK THE BAND
#############################################################

print(
  table(
    core_band$Phylum,
    useNA = "ifany"
  )
)


#############################################################
# REMOVE UNASSIGNED TIPS
#############################################################

core_band <- core_band[
  !is.na(core_band$Phylum),
  ,
  drop = FALSE
]


#############################################################
# ADD PHYLA BAND
#############################################################

p <- p +
  geom_fruit(
    data = core_band,
    geom = geom_tile,
    mapping = aes(
      y = label,
      x = 1,
      fill = Phylum
    ),
    offset = 0.02,
    pwidth = 0.045,
    width = 3
  )


#############################################################
# PHYLA COLORS
#############################################################

p <- p +
  scale_fill_manual(
    values = c(
      "Acidobacteriota"  = "#BCE1FF",
      "Actinobacteriota" = "#48C9B0",
      "Armatimonadota"   = "#A6A800",
      "Bacteroidota"     = "#EFB6D6",
      "Chloroflexi"      = "#919191",
      "Cyanobacteria"    = "#009051",
      "Firmicutes"       = "#F9FF91",
      "Gemmatimonadota"  = "#8C83E8",
      "Planctomycetota"  = "#F09263",
      "Proteobacteria"   = "#BDEB6E",
      "Verrucomicrobiota"   = "#941751",
      "Bdellovibrionota+Myxococcota" = "#73FDFF"
    ),
    na.value = "white"
  )

print(p)



ggsave(
  "LEfSe_genus_tree_2025_only.pdf",
  p,
  width = 10,
  height = 10,
  units = "in",
)


#############################################################
# EXPORT GENERA OVER-REPRESENTED IN BURN SITE TO EXCEL
#############################################################

library(writexl)

burn_genera_export <- results_filtered %>%
  dplyr::filter(Direction == "Burn") %>%
  dplyr::arrange(dplyr::desc(abs(LDA))) %>%
  dplyr::select(
    Genus,
    Family,
    Phylum,
    LDA,
    Direction,
    Burn_Mean,
    Control_Mean,
    Burn_Median,
    Control_Median
  )

writexl::write_xlsx(
  burn_genera_export,
  "Genera_Overrepresented_in_Burn.xlsx"
)

cat(
  "Exported",
  nrow(burn_genera_export),
  "Burn-overrepresented genera to:",
  file.path(getwd(), "Genera_Overrepresented_in_Burn.xlsx"),
  "\n"
)
