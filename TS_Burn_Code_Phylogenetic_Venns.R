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




#' @importFrom ggplot2 coord_cartesian
#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names prune_samples
#' @importFrom ape which.edge mrca nodepath
#' @importFrom utils combn
#' @export
#'
output_corePhyloVenn <- function(x,
                                 grouping,
                                 core_fraction = 0.5,
                                 ab_threshold1 = 0,
                                 ab_threshold2 = 0,
                                 ab_threshold3 = 0,
                                 mode = 'branch',
                                 selection = 'basic',
                                 max_tax=NULL,
                                 increase_cutoff = 2,
                                 initial_branches=NULL,
                                 rooted=TRUE,
                                 ordered_groups=NULL,
                                 show_percentage=TRUE,
                                 decimal=2,
                                 fill_color=c('red','orange','yellow','green','blue','purple','black'),
                                 fill_alpha=0.5,
                                 stroke_color='black',
                                 stroke_alpha = 1,
                                 stroke_size = 1,
                                 stroke_linetype = "solid",
                                 set_name_color = "black",
                                 set_name_size = 6,
                                 text_color = "black",
                                 text_size = 4,
                                 return_values = FALSE) {
  
  core<-core_fraction
  
  #Using an unrooted tree in branch mode with abundance thresholds is not currently supported.  If an unrooted tree is provided with non-zero abundance thresholds and branch mode is requested, reroot it
  if (rooted == FALSE){
    if (ab_threshold1+ab_threshold2+ab_threshold3>0){
      randoroot = sample(phy_tree(x)$tip.label, 1)
      warning("Use of an unrooted tree with non-zero abundance thresholds in branch mode is not supported.")
      warning("Randomly assigning root as -- ", randoroot, " -- in the phylogenetic tree provided.")
      phy_tree(x) <- root(phy=phy_tree(x), outgroup=randoroot, resolve.root=TRUE, interactive=FALSE)
      if( !is.rooted(phy_tree(x)) ){
        stop("Problem automatically rooting tree. Make sure your tree is rooted before attempting UniFrac calculation. See ?ape::root")
      }
      rooted = TRUE
    }
  }
  
  #find the number of different habitat types (e.g. hosts or environments) that are being compared
  group_count<-length(unique(grouping))
  
  #if no group order is specified, pick arbitrarily
  if (is.null(ordered_groups)){
    group_id<-unique(grouping)
    #otherwise use the specified group order
  }  else{group_id<-ordered_groups}
  
  #Venn diagrams can only be drawn for 7 or less habitats
  #If more than 7 habitats are entered, print out a warning
  if (group_count>7){
    warning('Warning: Too many habitat types!')
  }
  
  #otherwise proceed...
  else{
    
    #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
    newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)
    
    #if you are building a branch-based tree...
    if (mode=='branch'){
      
      #initialize a list of lists where the main list is of the habitats and the sublists are the samples associated with each habitat
      grouplist<-list()
      
      #initialize a list of lists where the main list is the habitats and the sublists are the edges associated with core microbes from each habitat
      edgelist<-list()
      
      #for each habitat...
      for (i in 1:group_count){
        #find the samples from that habitat
        temp<-list(sample_names(x)[which(grouping==group_id[i])])
        
        #put them into the list of samples in each habitat
        grouplist<-append(grouplist,temp)
        
        #Make a phyloseq object containing only samples from the focal habitat
        xg<-prune_samples(sample_names(x)[which(grouping==group_id[i])],x)
        
        #If using the basic selection mode...
        if (selection == 'basic'){
          
          #Call function to calculate all edges and core edges
          temp2<-basic_branch(xg,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)
          #If using the shade selection mode...
        }else if (selection == 'shade'){
          
          #If no maximum taxon is specified, default to considering 1%
          if (is.null(max_tax)){max_tax<-1}
          
          #Call the function to calculate all edges and core edges
          temp2<-shade_branch(xg,newtree,max_tax,increase_cutoff,initial_branches)
          
          #If a non-supported mode is selected print a warning and stop
        }else{
          
          stop('That method of core selection is not available. Please use basic or shade.')
          
        }
        
        #make a list of the core edges from the focal habitat
        core_edges<-temp2$core_edges
        
        #if the root is not included
        if (rooted==FALSE){
          
          #find the nodes associated with core edges
          nodes<-unique(c(newtree$edge[,1][core_edges],newtree$edge[,2][core_edges]))
          
          #find the mrca of each node in the tree
          cc<-mrca(newtree,full=TRUE)
          
          #find the mrca of each node associated with a core edge
          mrca_matrix<-cc[nodes,nodes]
          
          #find the unique mrcas for the core edge nodes
          mrca_list<-unique(as.vector(mrca_matrix))
          
          #find the unique mrcas plus core edge nodes
          mrca_list<-unique(mrca_list,nodes)
          
          #identify mrcas missing from the list of nodes associated with core edges
          missing<-mrca_list[which(!(mrca_list %in% nodes))]
          if (length(missing)>0){
            for (i in 1:length(missing)){
              for (j in 1:length(nodes)){
                #find the nodes connecting the missing mrcas to the nodes associated with core edges
                mrca_list<-c(mrca_list,nodepath(newtree,from=missing[i],to=nodes[j]))
              }
            }
          }
          
          #find the edges associated with all the nodes (core and mrcas)
          all_core_edges<-intersect(which(newtree$edge[,1] %in% mrca_list),which(newtree$edge[,2] %in% mrca_list))
          
          #add the missing edges to the core edges
          core_edges<-unique(c(core_edges,all_core_edges))
          
        }
        
        #append the edges from the focal habitat to the list of lists where the main list is of the habitats and the sublists are of the core edges for each habitat
        edgelist<-append(edgelist,list(core_edges))
        
      }
    }
    
    #if you are building a tip-based tree...
    else if (mode=='tip'){
      
      
      #initialize a list of lists where the main list is of the habitats and the sublists are the samples associated with each habitat
      grouplist<-list()
      
      #initialize a list of lists where the main list is the habitats and the sublists are the names of the core taxa associated with each habitat
      corelist<-list()
      
      #initialize a vector of all the names of the core taxa across any/all habitats
      allcorelist<-c()
      
      #for each habitat...
      for (i in 1:group_count){
        
        #find the samples from that habitat
        temp<-list(sample_names(x)[which(grouping==group_id[i])])
        
        #put them into the list of samples in the habitat list
        grouplist<-append(grouplist,temp)
        
        #make a phyloseq object for habitat
        xg<-prune_samples(sample_names(x)[which(grouping==group_id[i])],x)
        
        if (selection == 'basic'){
          
          #Call function to calculate all edges and core edges
          temp2<-basic_tip(xg,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)
          
        }else if (selection == 'shade'){
          
          #If no maximum taxon is specified, default to considering 1%
          if (is.null(max_tax)){max_tax<-1}
          temp2<-shade_tip(xg,newtree,max_tax,increase_cutoff)
          
        }else{
          
          warning('That method of core selection is not available. Please use basic or shade.')
          
        }
        
        #put them into the list of taxa in the habitat list
        corelist<-append(corelist,list(temp2$core_taxa))
        
        #put them into the vector of core taxa from any/all habitats
        allcorelist<-unique(c(temp$core_taxa,allcorelist))
        
      }
      
      #find the edges associated with core taxa from any/all habitats (this is the minimal spanning tree)
      spanlist<-which.edge(newtree,allcorelist)
      habitatspanlist<-c()
      for (i in 1:group_count){
        habitatspanlist<-c(habitatspanlist,which.edge(newtree,corelist[[i]]))
      }
      
      #make a list of lists, with the main list being habitats and the sublists being the edges associated with core taxa in each habitat
      edgelist<-list()
      for (i in 1:group_count){
        #initially, include the edges associated with the root
        edgelist<-append(edgelist,list(which.edge(newtree,c('outgroup',corelist[[i]]))))
      }
      
      
      #if you are including the root...
      if (rooted ==TRUE){
        #if you are not including the root...
      }else{
        #for each habitat...
        for (i in 1:group_count){
          #remove all edges that are not part of the minimal spanning tree for core microbes from each habitat
          edgelist[[i]]<-edgelist[[i]][which(edgelist[[i]] %in% habitatspanlist[[i]])]
        }
      }
      
      #if a mode that is not supported is entered, print a warning
    }else{stop('Warning: that mode is not supported')}
    
    #initialize a list of possible habitat combinations
    combos<-list()
    
    #initialize a list of lists, with the main lists being habitat combinations and the sublists being the branch lengths shared by the habitat combinations
    intersections<-list()
    
    #initialize a list of total branch length shared by each habitat combination
    lengths<-c()
    
    #for anywhere from 1 to n combinations (i.e., {1,2} is a length 2 combination, {1,4,5} is a length 3 combination), where n is the total number of habitats...
    for (i in 1:group_count){
      
      
      #list all of the possible combinations of that size (e.g., {1,2},{1,3},{2,3})
      ff<-combn(group_count,i)
      
      #count the possible combinations of that size
      no_combinations<-length(ff[1,])
      
      #find the length of the combination
      combination_length<-length(ff[,1])
      
      #for each possible combination of the focal size...
      for (j in 1:no_combinations){
        
        #add the particular combination of habitats to your list of possible habitat combinations
        combos<-append(combos,list(ff[,j]))
        
        #find the habitats that are outside the particular combination
        outside<-which(!(1:1:group_count %in% ff[,j]))
        
        #find the edges in the first habitat of the combination
        shared_temp<-edgelist[[ff[1,j]]]
        
        
        #if there are more habitats in the combination...
        if (combination_length>1){
          
          #find the edges that are shared by all habitats in the combination
          for (k in 2:combination_length){
            shared_temp<-intersect(shared_temp,edgelist[[ff[k,j]]])
          }
          
        }
        
        #if there are any habitats outside the combination...
        if (length(outside)>0){
          
          #for each habitat outside...
          for (g in 1:length(outside)){
            #remove the edges that it shares with the focal combination of habitats (shared edges must be exclusive to the focal habitat combination)
            shared_temp<-setdiff(shared_temp,edgelist[[outside[g]]])
          }
          
        }
        
        #add the edges that are exclusively shared with the focal combination of habitats to the list of shared edges for each habitat combination
        intersections<-append(intersections,list(shared_temp))
        
        
        #add the lengths of the branches shared with the focal combination of habitats to the list of lengths for each habitat combination
        if (rooted==TRUE){
          
          #if the root was included, then you need to remove the branch length for the added outgroup that was shared with everyone
          if (length(shared_temp)>0){
            lengths<-c(lengths,sum(newtree$edge.length[shared_temp])-0.0001)
          }else{lengths<-c(lengths,0)}
          
        }else{
          
          #if the root was not included, then there is no correction to branch lengths
          if (length(shared_temp)>0){
            lengths<-c(lengths,sum(newtree$edge.length[shared_temp]))
          }else{lengths<-c(lengths,0)}
          
        }
      }
    }
    
    if (sum(lengths)==0){stop('There are no core taxa in any group')}
    
    # Output values for each exclusive Venn compartment
    compartment_values <- data.frame(
      Compartment = vapply(
        combos,
        function(idx) paste(group_id[idx], collapse = " & "),
        character(1)
      ),
      Branch_length = lengths,
      Percentage = 100 * lengths / sum(lengths)
    )
    
    compartment_values$Branch_length <-
      round(compartment_values$Branch_length, decimal)
    
    compartment_values$Percentage <-
      round(compartment_values$Percentage, decimal)
    
    
    
    #if the Venn diagram is shown as a percentage of the branch lengths...
    if (show_percentage){
      #round to the decimal desired, multiple by 100 and divide by the total length
      fractions<-round(10^(decimal+2)*lengths/sum(lengths))
    }else{
      #otherwise just round to the decimal desired
      fractions<-round(10^(decimal)*lengths)
    }
    #note that the above are in whole numbers... so if you want 25.56% you have 2556 in that slot...
    #this is for making the dummy dataframe that can be read by the Venn diagram package
    #find the cumulative percentages or branch lengths across the list of possible habitat combinations
    cumfractions<-cumsum(fractions)
    
    #initialize a dataframe for the Venn diagram package with FALSE for every entry
    #This dataframe has as many rows as the total value of the cumulative percentages or branch lengths in your cumfractions vector
    #So, for instance, if you want percentages with 0 decimals, you'd have 100 rows of FALSE entries... if you want percentages with 2 decimals you'd have 10000 rows of FALSE entries
    final<-data.frame(rep(FALSE,cumfractions[length(cumfractions)]))
    #There should be as many columns as there are habitats
    for (i in 2:group_count){
      final<-cbind(final,rep(FALSE,cumfractions[length(cumfractions)]))
    }
    #for each row in the dataframe
    counter<-1
    for (i in 1:length(cumfractions)){
      if (counter<cumfractions[i]){
        #find the habitat combinations associated with those percentages and set them to TRUE
        final[counter:cumfractions[i],combos[i][[1]]]<-TRUE
        counter<-cumfractions[i]+1
      }
    }
    
    #name the dataframe columns based on the habitats
    colnames(final)<-group_id
    #Plot the Venn diagram
    venn_plot <- ggvenn2(
      final, columns = NULL, show_elements = FALSE,
      show_percentage, digits = decimal, fill_color, fill_alpha,
      stroke_color, stroke_alpha, stroke_size, stroke_linetype,
      set_name_color, set_name_size, text_color, text_size
    ) + coord_cartesian(clip = "off")
    
    attr(venn_plot, "compartment_values") <- compartment_values
    
    return(venn_plot)  }
}

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
Z_physeq <- prune_samples(sample_data(Z_physeq)$QIIME_ID_aligned != 'C1-25-27-NoBurn', Z_physeq)
Z_physeq <- prune_samples(sample_data(Z_physeq)$QIIME_ID_aligned != 'C1-25-27-NoBurn', Z_physeq)
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

#Rename Sites
sample_data(rZOTU)$Site[which(sample_data(rZOTU)$Site == '31_burn')]<-'Burn'
sample_data(rZOTU)$Site[which(sample_data(rZOTU)$Site == '32_no_burn')]<-'Control_2'
sample_data(rZOTU)$Site[which(sample_data(rZOTU)$Site == '2950')]<-'Control_1'

SitexYear<-c()
Year<-c()
for (k in 1:length(sample_data(rZOTU)$Site)){
  SitexYear<-c(SitexYear,paste0(sample_data(rZOTU)$Site[k],'_',strsplit(sample_data(rZOTU)$Season[k],'_')[[1]][2])) 
  Year<-c(Year,strsplit(sample_data(rZOTU)$Season[k],'_')[[1]][2])
}

sample_data(rZOTU)$SitexYear<-SitexYear
sample_data(rZOTU)$Year<-Year

sample_data(rZOTU)$BurnStatus <- ifelse(
  sample_data(rZOTU)$Site == "Burn",
  "Burn",
  "Control"
)

sample_data(rZOTU)$BurnStatus <- factor(
  sample_data(rZOTU)$BurnStatus,
  levels = c("Control", "Burn")
)


n_trials <- 10
core_fraction <- 0.9

############################
# SETTINGS
############################

set.seed(12345)

n_per_site <- 14

ordered_groups <- c("Burn", "Control_2", "Control_1")
fill_colors <- c("red", "green", "forestgreen")


############################
# PREPARE 2025 DATA
############################

rZOTU_2025 <- prune_samples(
  sample_data(rZOTU)$Season != "Spring_2024",
  rZOTU
)

# Split samples by site
samples_by_site <- split(
  sample_names(rZOTU_2025),
  sample_data(rZOTU_2025)$Site
)

# Confirm that each site has enough samples
print(lengths(samples_by_site))

if (any(lengths(samples_by_site[ordered_groups]) < n_per_site)) {
  stop("At least one site has fewer than 14 available samples.")
}


############################
# RUN 10 TRIALS
############################

trial_values <- vector("list", n_trials)

for (trial in seq_len(n_trials)) {
  
  message("Running trial ", trial, " of ", n_trials)
  
  # Randomly sample 14 samples from each site
  keep <- unlist(
    lapply(
      samples_by_site[ordered_groups],
      function(samples) sample(samples, n_per_site)
    ),
    use.names = FALSE
  )
  
  rZOTU_trial <- prune_samples(keep, rZOTU_2025)
  
  # Generate the phylogenetic Venn diagram
  trial_plot <- output_corePhyloVenn(
    rZOTU_trial,
    sample_data(rZOTU_trial)$Site,
    core_fraction,
    ordered_groups = ordered_groups,
    fill_color = fill_colors,
    show_percentage = FALSE
  )
  
  # Extract compartment values
  values <- attr(trial_plot, "compartment_values")
  
  if (is.null(values)) {
    stop(
      "Trial ", trial,
      " did not return the 'compartment_values' attribute."
    )
  }
  
  # Save results with trial number
  values$Trial <- trial
  trial_values[[trial]] <- values
}


############################
# COMBINE TRIAL RESULTS
############################

all_trial_values <- do.call(rbind, trial_values)

# Average branch length for each compartment
average_values <- aggregate(
  Branch_length ~ Compartment,
  data = all_trial_values,
  FUN = mean
)

# Calculate standard deviation across trials
sd_values <- aggregate(
  Branch_length ~ Compartment,
  data = all_trial_values,
  FUN = sd
)

names(sd_values)[2] <- "SD"

# Combine mean and standard deviation
average_values <- merge(
  average_values,
  sd_values,
  by = "Compartment",
  sort = FALSE
)

# Add percentage of the mean total branch length
average_values$Percentage <-
  100 * average_values$Branch_length /
  sum(average_values$Branch_length)

average_values$Branch_length <-
  round(average_values$Branch_length, 6)

average_values$SD <- round(average_values$SD, 6)
average_values$Percentage <- round(average_values$Percentage, 2)

print(average_values)

############################
# PLOT AVERAGE VENN DIAGRAM
############################

# Use the compartment order expected by the original function
compartments <- c(
  "Burn",
  "Control_2",
  "Control_1",
  "Burn & Control_2",
  "Burn & Control_1",
  "Control_2 & Control_1",
  "Burn & Control_2 & Control_1"
)

# Match averaged values to the expected compartment order
mean_lengths <- setNames(
  average_values$Branch_length,
  average_values$Compartment
)

    
# Make sure the compartment values are in the correct order
mean_lengths <- mean_lengths[compartments]

scale_factor <- 100
# Convert to integer row counts
plot_counts <- round(mean_lengths * scale_factor)

# Create exactly enough rows for all compartments
final_avg <- data.frame(
  Burn = rep(FALSE, sum(plot_counts)),
  Control_2 = rep(FALSE, sum(plot_counts)),
  Control_1 = rep(FALSE, sum(plot_counts))
)

# Assign every row to its exclusive compartment
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

# Check: every row should belong to at least one group
stopifnot(all(rowSums(final_avg) > 0))

# Plot
p_average <- ggvenn2(
  final_avg,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = FALSE,
  digits = 2,
  fill_color = c("red", "green", "forestgreen"),
  fill_alpha = 0.5,
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 6,
  text_color = "black",
  text_size = 8
) + ggplot2::coord_cartesian(clip = "off")

print(p_average)

############################
# SAVE RESULTS
############################

write.csv(
  all_trial_values,
  "2025_phylo_venn_10_trials.csv",
  row.names = FALSE
)

write.csv(
  average_values,
  "2025_phylo_venn_averages.csv",
  row.names = FALSE
)




############################
# SETTINGS
############################

set.seed(12345)

n_per_site <- 10


############################
# PREPARE 2024 DATA
############################

rZOTU_2025 <- prune_samples(
  sample_data(rZOTU)$Season != "Spring_2024",
  rZOTU
)

rZOTU_2024 <- prune_samples(
  sample_data(rZOTU)$Season != "Spring_2025",
  rZOTU
)

# Split samples by site
samples_by_site <- split(
  sample_names(rZOTU_2024),
  sample_data(rZOTU_2024)$Site
)

# Confirm sample counts
print(lengths(samples_by_site))

if (any(lengths(samples_by_site[ordered_groups]) < n_per_site)) {
  stop("At least one site has fewer than 10 available samples.")
}


############################
# RUN 10 TRIALS
############################

trial_values <- vector("list", n_trials)

for (trial in seq_len(n_trials)) {
  
  message("Running trial ", trial, " of ", n_trials)
  
  # Randomly sample 10 samples from each site
  keep <- unlist(
    lapply(
      samples_by_site[ordered_groups],
      function(samples) sample(samples, n_per_site)
    ),
    use.names = FALSE
  )
  
  rZOTU_trial <- prune_samples(keep, rZOTU_2024)
  
  # Generate the phylogenetic Venn diagram
  trial_plot <- output_corePhyloVenn(
    rZOTU_trial,
    sample_data(rZOTU_trial)$Site,
    core_fraction,
    ordered_groups = ordered_groups,
    fill_color = fill_colors,
    show_percentage = FALSE
  )
  
  # Extract compartment values
  values <- attr(trial_plot, "compartment_values")
  
  if (is.null(values)) {
    stop(
      "Trial ", trial,
      " did not return the 'compartment_values' attribute."
    )
  }
  
  values$Trial <- trial
  trial_values[[trial]] <- values
}


############################
# COMBINE TRIAL RESULTS
############################

all_trial_values <- do.call(rbind, trial_values)

# Mean branch length for each compartment
average_values <- aggregate(
  Branch_length ~ Compartment,
  data = all_trial_values,
  FUN = mean
)

# Standard deviation across trials
sd_values <- aggregate(
  Branch_length ~ Compartment,
  data = all_trial_values,
  FUN = sd
)

names(sd_values)[2] <- "SD"

average_values <- merge(
  average_values,
  sd_values,
  by = "Compartment",
  sort = FALSE
)

# Percentage of total mean branch length
average_values$Percentage <-
  100 * average_values$Branch_length /
  sum(average_values$Branch_length)

# Round values
average_values$Branch_length <-
  round(average_values$Branch_length, 6)

average_values$SD <- round(average_values$SD, 6)
average_values$Percentage <- round(average_values$Percentage, 2)


############################
# ORDER AND DISPLAY RESULTS
############################

compartments <- c(
  "Burn",
  "Control_2",
  "Control_1",
  "Burn & Control_2",
  "Burn & Control_1",
  "Control_2 & Control_1",
  "Burn & Control_2 & Control_1"
)

average_values$Compartment <- factor(
  average_values$Compartment,
  levels = compartments
)

average_values <- average_values[
  order(average_values$Compartment),
]

print(average_values)


############################
# SAVE RESULTS
############################

write.csv(
  all_trial_values,
  "2024_phylo_venn_10_trials.csv",
  row.names = FALSE
)

write.csv(
  average_values,
  "2024_phylo_venn_averages.csv",
  row.names = FALSE
)

############################
# PLOT AVERAGE 2024 VENN DIAGRAM
############################

# Define compartment order
compartments <- c(
  "Burn",
  "Control_2",
  "Control_1",
  "Burn & Control_2",
  "Burn & Control_1",
  "Control_2 & Control_1",
  "Burn & Control_2 & Control_1"
)

# Extract mean branch lengths in the correct order
mean_lengths <- setNames(
  average_values$Branch_length,
  as.character(average_values$Compartment)
)

mean_lengths <- mean_lengths[compartments]
mean_lengths[is.na(mean_lengths)] <- 0

# Use the same scaling approach as your working 2025 plot
scale_factor <- 100
plot_counts <- round(mean_lengths * scale_factor)

# Create dummy data frame with exactly enough rows
final_avg <- data.frame(
  Burn = rep(FALSE, sum(plot_counts)),
  Control_2 = rep(FALSE, sum(plot_counts)),
  Control_1 = rep(FALSE, sum(plot_counts))
)

# Assign each row to its exclusive compartment
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

# Confirm no rows are outside all groups
stopifnot(all(rowSums(final_avg) > 0))

# Plot
p_average_2024 <- ggvenn2(
  final_avg,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = FALSE,
  digits = 2,
  fill_color = c("red", "green", "forestgreen"),
  fill_alpha = 0.5,
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 6,
  text_color = "black",
  text_size = 8
) + ggplot2::coord_cartesian(clip = "off")

print(p_average_2024)


############################
# SETTINGS
############################

set.seed(12345)



############################
# PREPARE 2024 AND 2025 DATA
############################


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


############################
# RUN 10 TRIALS
############################

trial_values <- vector("list", n_trials)

for (trial in seq_len(n_trials)) {
  
  message("Running trial ", trial, " of ", n_trials)
  
  # Sample 14 from each site in 2025
  keep1 <- unlist(
    lapply(
      samples_2025_by_site[ordered_groups],
      function(x) sample(x, 14)
    ),
    use.names = FALSE
  )
  
  # Sample 10 from each site in 2024
  keep2 <- unlist(
    lapply(
      samples_2024_by_site[ordered_groups],
      function(x) sample(x, 10)
    ),
    use.names = FALSE
  )
  
  # Combine samples from both years
  keep <- c(keep1, keep2)
  
  rZOTU_trial <- prune_samples(keep, rZOTU)
  
  # Generate phylogenetic Venn diagram
  trial_plot <- output_corePhyloVenn(
    rZOTU_trial,
    sample_data(rZOTU_trial)$Site,
    core_fraction,
    ordered_groups = ordered_groups,
    fill_color = fill_colors,
    show_percentage = FALSE
  )
  
  # Extract compartment values
  values <- attr(trial_plot, "compartment_values")
  
  if (is.null(values)) {
    stop(
      "Trial ", trial,
      " did not return the 'compartment_values' attribute."
    )
  }
  
  values$Trial <- trial
  trial_values[[trial]] <- values
}


############################
# COMBINE TRIAL RESULTS
############################

all_trial_values <- do.call(rbind, trial_values)

# Average branch length for each compartment
average_values <- aggregate(
  Branch_length ~ Compartment,
  data = all_trial_values,
  FUN = mean
)

# Standard deviation across trials
sd_values <- aggregate(
  Branch_length ~ Compartment,
  data = all_trial_values,
  FUN = sd
)

names(sd_values)[2] <- "SD"

average_values <- merge(
  average_values,
  sd_values,
  by = "Compartment",
  sort = FALSE
)

# Percentage of total mean branch length
average_values$Percentage <-
  100 * average_values$Branch_length /
  sum(average_values$Branch_length)

average_values$Branch_length <-
  round(average_values$Branch_length, 6)

average_values$SD <- round(average_values$SD, 6)
average_values$Percentage <- round(average_values$Percentage, 2)


############################
# ORDER AND DISPLAY RESULTS
############################

compartments <- c(
  "Burn",
  "Control_2",
  "Control_1",
  "Burn & Control_2",
  "Burn & Control_1",
  "Control_2 & Control_1",
  "Burn & Control_2 & Control_1"
)

average_values$Compartment <- factor(
  average_values$Compartment,
  levels = compartments
)

average_values <- average_values[
  order(average_values$Compartment),
]

print(average_values)


############################
# SAVE RESULTS
############################

write.csv(
  all_trial_values,
  "2024_2025_phylo_venn_10_trials.csv",
  row.names = FALSE
)

write.csv(
  average_values,
  "2024_2025_phylo_venn_averages.csv",
  row.names = FALSE
)


############################
# PLOT AVERAGE VENN DIAGRAM
############################

mean_lengths <- setNames(
  average_values$Branch_length,
  as.character(average_values$Compartment)
)

mean_lengths <- mean_lengths[compartments]
mean_lengths[is.na(mean_lengths)] <- 0

# Convert branch lengths to plotting counts,
# using the same scaling as your working plots
scale_factor <- 100
plot_counts <- round(mean_lengths * scale_factor)

# Create dummy data frame with exactly enough rows
final_avg <- data.frame(
  Burn = rep(FALSE, sum(plot_counts)),
  Control_2 = rep(FALSE, sum(plot_counts)),
  Control_1 = rep(FALSE, sum(plot_counts))
)

# Assign every row to its exclusive compartment
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

# Confirm no empty compartments are accidentally added
stopifnot(all(rowSums(final_avg) > 0))

# Plot
p_average_2024_2025 <- ggvenn2(
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
  set_name_size = 10,
  text_color = "black",
  text_size = 8
) + ggplot2::coord_cartesian(clip = "off")

print(p_average_2024_2025)





