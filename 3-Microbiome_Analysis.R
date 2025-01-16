# MADE BY BRYAN LE, NOVEMBER 2024

# This R script builds on the output files from the pipeline to perform in-depth
# statistical analyses. It examines taxonomic differences, alpha diversity 
# differences, beta diversity differences, and functional differences between 
# healthy individuals and those with inflammatory bowel disease (IBD). The goal 
# is to identify and characterize variations in the core microbiomes of the two 
# groups, providing insights into the microbiological distinctions associated 
# with IBD.

# The taxon level is currently set as 'L7' (species).  To perform analysis at a 
# different, find all instances of 'L7' and replace it with the desired level. 



################################################################################
### PART 0: LOADING MODULES AND DATA, AND DEFINING PATHWAYS
library('vegan')
library('car')
library('ggplot2')
library('tidyr')
library('Maaslin2')

OUTPUT_DIR <- "INSERT PATHWAY TO OUTPUT DIRECTORY HERE"

meta_data <- read.delim(paste0(OUTPUT_DIR, "/meta_data/metadata.txt"),row=1, as.is=FALSE)
alpha <- read.delim(paste0(OUTPUT_DIR, "/Kraken/alpha/L7_alpha-diversity.txt"),row=1)
beta_bj <- read.delim(paste0(OUTPUT_DIR, "/Kraken/beta/L7_beta/binary_jaccard_taxa_table_L7_final.txt"),row=1)
species <- read.delim(paste0(OUTPUT_DIR, "/Kraken/Taxa/taxa_table_L7_final_transposed_relative_abundance.txt"),row=1)
paths <- read.delim(paste0(OUTPUT_DIR, "/HUMAaN/Merged_Files/normalized_pathabundance.tsv"), row=1)

meta_data$Condition <- factor(meta_data$Condition,levels=c('Healthy','IBD'))
beta_bj <- as.matrix(beta_bj)
map <- meta_data[row.names(species), , drop = FALSE]



################################################################################
######### PART 1: Exploratory Analysis - Going through all species and using ANOVA to test for relative OTU abundance association with health conditions ### 

otu_results <- data.frame(Species = character(), p_value = numeric(), adjusted_p_value = numeric(), stringsAsFactors = FALSE)
for (i in seq_along(colnames(species))) {
  species_name <- colnames(species)[i]
  
  model <- lm(species[[species_name]] ~ map$Condition)
  anova_result <- anova(model)
  
  p_value <- anova_result$`Pr(>F)`[1]
  
  otu_results <- rbind(otu_results, data.frame(Species = species_name, p_value = p_value, stringsAsFactors = FALSE))
}

otu_results$adjusted_p_value <- p.adjust(otu_results$p_value, method = "fdr")
significant_otus <- otu_results[otu_results$p_value < 0.05, ]

# This is how many OTUS are significantly associated with our groups before FDR correction
print(length(significant_otus$Species))



################################################################################
### PART 2: Visualizing differences in average relative abundance of top 20 significant OTUS between health and IBD individuals ### 
top_20_otus <- significant_otus[order(significant_otus$p_value)[1:20], ]
subset_species <- species[, top_20_otus$Species]
subset_species_healthy <- subset_species[rownames(map)[map$Condition == "Healthy"], ]
subset_species_IBD <- subset_species[rownames(map)[map$Condition == "IBD"], ]

average_abundance <- data.frame(Species = character(), average_relative_abundance = numeric(), condition = character(), stringsAsFactors = FALSE)
for (i in top_20_otus$Species) {
  average_relative_abundance_healthy <- mean(subset_species_healthy[, i], na.rm = TRUE)
  average_relative_abundance_ibd <- mean(subset_species_IBD[, i], na.rm = TRUE)
  
  new_row_healthy <- data.frame(
    Species = i, 
    average_relative_abundance = average_relative_abundance_healthy, 
    condition = "Healthy", 
    stringsAsFactors = FALSE
  )
  
  new_row_ibd <- data.frame(
    Species = i, 
    average_relative_abundance = average_relative_abundance_ibd, 
    condition = "IBD", 
    stringsAsFactors = FALSE
  )
    average_abundance <- rbind(average_abundance, new_row_healthy)
    average_abundance <- rbind(average_abundance, new_row_ibd)
}


set.seed(333333)
random_colors <- sample(colors(), length(unique(average_abundance$Species)))
ggplot(average_abundance, aes(x = condition, y = average_relative_abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(x = "Condition", y = "Average Relative Abundance") +
  theme_classic() +
  scale_fill_manual(values = random_colors) +
  theme(legend.key.size = unit(0.1, "cm"), legend.text = element_text(size = 8), 
    legend.title = element_text(size = 9) 
  )



################################################################################
### PART 3: Visualizing differences in OTU association before and after FDR p-value correction
otu_results$p_value_rank <- rank(otu_results$p_value, ties.method = "first")
otu_results$corrected_p_value_rank <- rank(otu_results$adjusted_p_value, ties.method = "first")
otu_results_long <- otu_results %>%
  pivot_longer(cols = c(p_value, adjusted_p_value), names_to = "Value_Type",values_to = "Value"
  )

GROUP.COLORS.FADED <- c("p_value" = "#00F494", "adjusted_p_value" = "#C77CFF", "Significance Cutoff" = "black")
ggplot(otu_results_long, aes(x = p_value_rank, y = Value, color = Value_Type)) +
  geom_line(size = 2) +
  geom_point(size = 0) +
  geom_line(aes(y = 0.05, color = "Significance Cutoff"), linetype = "dashed", size = 1) +
  scale_color_manual(values = GROUP.COLORS.FADED) +
  scale_linetype_manual(values = c("Significance Cutoff" = "dashed")) +
  labs(x = "Rank", y = "P-value") +
  theme_classic()



################################################################################
### PART 4: Exploratory Analysis - Finding Core Microbiome of Healthy and IBD Individuals; Finding OTUS Exclusive to Either Group

# Finding OTUs Present in all Healthy individuals
otus_in_healthy_individuals <- species[map$Condition == "Healthy", , drop = FALSE]
otus_present_in_all_healthy_individuals <- otus_in_healthy_individuals[, 
                                                                       which(apply(otus_in_healthy_individuals, 2, function(col) all(col > 0))), 
                                                                       drop = FALSE
]
healthy_otus <- colnames(otus_present_in_all_healthy_individuals)


# Finding OTUs Present in all IBD individuals
otus_in_ibd_individuals <- species[map$Condition == "IBD", , drop = FALSE]
otus_present_in_all_ibd_individuals <- otus_in_ibd_individuals[, 
                                                               which(apply(otus_in_ibd_individuals, 2, function(col) all(col > 0))), 
                                                               drop = FALSE
]
ibd_otus <- colnames(otus_present_in_all_ibd_individuals)


# Finding OTUs Exclusive to IBD and Healthy individuals
ibd_exclusive_otus <- setdiff(ibd_otus, colnames(otus_in_healthy_individuals))
healthy_exclusive_otus <- setdiff(healthy_otus, colnames(otus_in_ibd_individuals))

# This is the number of OTUS present in all healthy individuals
print(length(healthy_otus))
# This is the number of OTUS present in all IBD individuals
print(length(ibd_otus))
# This is the number of OTUs exclusive to healthy individuals
print(length(healthy_exclusive_otus))
# This is the number of OTUs exclusive to IBD individuals
print(length(ibd_exclusive_otus))

healthy_averages <- apply(otus_present_in_all_healthy_individuals, 2, mean, na.rm = TRUE)
healthy_averages_df <- data.frame(Column = names(healthy_averages), relative_abundance = healthy_averages)
healthy_averages_subset <- healthy_averages_df[order(-healthy_averages_df$relative_abundance)[1:10], ]

ibd_averages <- apply(otus_present_in_all_ibd_individuals, 2, mean, na.rm = TRUE)
ibd_averages_df <- data.frame(Column = names(ibd_averages), relative_abundance = ibd_averages)
ibd_averages_subset <- ibd_averages_df[order(-ibd_averages_df$relative_abundance)[1:10], ]

ggplot(healthy_averages_subset, aes(x = reorder(Column, -relative_abundance), y = relative_abundance, fill = Column)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "A", x = "Species", y = "Average Relative Abundance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.95, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggplot(ibd_averages_subset, aes(x = reorder(Column, -relative_abundance), y = relative_abundance, fill = Column)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "B", x = "Species", y = "Average Relative Abundance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.95, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
  )


ggplot(ibd_diff_subset, aes(x = reorder(Species, -diff_value), y = diff_value, fill = Species)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "B", x = "Species", y = "Difference in Average Relative Abundance") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.95, face = "bold"),  
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
  )



################################################################################
### PART 5: Exploratory Analysis - Finding OTUs that are most differentially abundant in both groups
subset_species_healthy <- species[rownames(map)[map$Condition == "Healthy"], ]
subset_species_IBD <- species[rownames(map)[map$Condition == "IBD"], ]
diff_average_abundance <- data.frame(Species = character(), diff_value = numeric(), condition = character(), stringsAsFactors = FALSE)
for (i in colnames(species)) {
  average_relative_abundance_healthy <- mean(subset_species_healthy[, i], na.rm = TRUE)
  average_relative_abundance_ibd <- mean(subset_species_IBD[, i], na.rm = TRUE)
  
  new_row_healthy <- data.frame(
    Species = i, 
    diff_value = average_relative_abundance_healthy - average_relative_abundance_ibd, 
    condition = "Healthy", 
    stringsAsFactors = FALSE
  )
  
  new_row_ibd <- data.frame(
    Species = i, 
    diff_value = average_relative_abundance_ibd - average_relative_abundance_healthy, 
    condition = "IBD", 
    stringsAsFactors = FALSE
  )
  diff_average_abundance <- rbind(diff_average_abundance, new_row_healthy)
  diff_average_abundance <- rbind(diff_average_abundance, new_row_ibd)
}

healthy_diff <- diff_average_abundance[diff_average_abundance$condition=="Healthy", ]
ibd_diff <- diff_average_abundance[diff_average_abundance$condition=="IBD", ]

healthy_diff_subset <- healthy_diff[order(-healthy_diff$diff_value)[1:10], ]
ibd_diff_subset <- ibd_diff[order(-ibd_diff$diff_value)[1:10], ]

ggplot(healthy_diff_subset, aes(x = reorder(Species, -diff_value), y = diff_value, fill = Species)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "A", x = "Species", y = "Difference in Average Relative Abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.95, face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


ggplot(ibd_diff_subset, aes(x = reorder(Species, -diff_value), y = diff_value, fill = Species)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "B", x = "Species", y = "Difference in Average Relative Abundance") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.95, face = "bold"),  
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



################################################################################
#### PART 6: BETA ANALYSIS
pc <- cmdscale(as.dist(beta_bj))
pcoa_df <- data.frame(
  Sample = rownames(pc),
  PC1 = pc[, 1],
  PC2 = pc[, 2],
  Group = map$Condition
)

GROUP.COLORS.FADED <- c("IBD" = "coral1", "Healthy" = "deepskyblue")
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  xlab("PC1") +
  ylab("PC2") +
  scale_color_manual(values = GROUP.COLORS.FADED) +
  stat_ellipse(type = "t", level = 0.68, size = 1,alpha = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    axis.line = element_line(color = "black", linewidth = 0.005)
  )

# PERMANOVA association testing for Binary-Jaccard beta diversity
adonis2(beta_bj ~ Condition, data=map)



################################################################################
### PART 7: FUNCTIONAL ANALYSIS

# Running MaAsLin2 to find pathways enriched in IBD individuals (relative to healthy individuals)
fit <- Maaslin2(
  input_data = paths,      
  input_metadata = map,       
  output = paste0(OUTPUT_DIR, "/Maaslin2/Results"),    
  fixed_effects = c("HEALTH_STATUS"),      
  random_effects = NULL,           
  normalization = "NONE",          
  transform = "LOG"              
)

all_results <- read.table(paste0(OUTPUT_DIR, "/Maaslin2/Results"), header = TRUE, sep = "\t")

# Cleaning up data
top_results <- all_results[order(all_results$qval), ]
subset_results <- top_results[, grep("feature|value|pval|qval", colnames(top_results), ignore.case = TRUE)]
