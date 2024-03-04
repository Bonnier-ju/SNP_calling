# Define the base path for input files
setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.2-normal_filtering")
base_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.2-normal_filtering"

# Path where saving plot
plot_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.2-normal_filtering/Plots/"


# Define the path for the population information CSV file
population_info_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.2-normal_filtering/Population_info.csv"


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)
library(data.table)
library(scales) 
library(ggrepel)
library(tidyverse)



######################### PCA analysis for individuals #############################
####################################################################################

# Read PCA results and population information
pca_data <- read_delim(file.path(base_path, "missing_filtered_stats.eigenvec"), 
                       delim = "\t") # Read PCA results with a tab delimiter
population_info <- read_csv(population_info_path) # Read population info, assuming CSV format

#Variance %
eigenvalues <- read.table("missing_filtered_stats.eigenval", header = FALSE)
eigenvalues <- eigenvalues$V1

total_variance <- sum(eigenvalues)
percent_variance <- (eigenvalues / total_variance) * 100

pc1_percent <- percent_variance[1]
pc2_percent <- percent_variance[2]
pc3_percent <- percent_variance[3]

# Prepare PCA data for merging by renaming IID to SampleID
# and creating a new column with the first 7 letters of each SampleID
pca_data <- pca_data %>%
  rename(SampleID = IID) %>%
  mutate(SampleID_short = substr(SampleID, 1, 7)) # Truncate SampleID to first 7 letters

# Merge PCA data with population information based on SampleID
pca_data <- pca_data %>%
  inner_join(population_info, by = "SampleID")

pop_order <- c("Apatou", "Acarouany", "Piste_St_Elie", "Nouragues_Inselberg", 
               "Montagne_tortue", "Cacao", "Route_de_l'est", "Foret_Regina_St_Georges", "Regina", "Saut_Lavilette", "St_georges")

pca_data$Pop <- factor(pca_data$Pop, levels = pop_order)


my_colors <- c("Apatou" = "#8B2500", "Acarouany" = "#CD3700",  "Piste_St_Elie" = "#FF4500",
               "Nouragues_Inselberg" = "darkolivegreen1", "Montagne_tortue" = "darkolivegreen3","Cacao" = "darkslategray2",
               "Route_de_l'est" = "#00BFFF", "Foret_Regina_St_Georges" = "deepskyblue3","Regina" = "#1C86EE",
               "Saut_Lavilette" = "#1874CD",  "St_georges" = "dodgerblue4")

# Visualisation PCA avec ggplot2, en utilisant la palette de couleurs personnalisée
ggplot(pca_data, aes(x = PC2, y = PC3, label = SampleID_short, color = Pop)) +
  geom_point() + 
  geom_text_repel(size = 3) +
  #stat_ellipse(type = "t", level = 0.95) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "PCA Plot", 
       x = paste("PC2 (", round(pc2_percent, 2), "%)", sep=""), 
       y = paste("PC3 (", round(pc3_percent, 2), "%)", sep=""),
       color = "Population") + 
  scale_color_manual(values = my_colors) + # Utiliser les couleurs définies manuellement
  theme(legend.position = "right", legend.title = element_text(face = "bold"))



############################################## Heterozygosity rate per sample ######################################################
####################################################################################################################################

het_data <- read.table("missing_filtered_stats.het", header = TRUE)

# Truncate individual IDs to the first 7 characters
het_data$IID_short <- substr(het_data$IID, 1, 7)


#O(HOM) (Observed Homozygotes): The observed number of homozygous genotypes in the individual. 
#It is the total number of loci where the individual is homozygous.

ggplot(het_data, aes(x = IID_short, y = `O.HOM.`)) + 
  geom_bar(stat = "identity", fill = "deeppink4") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Heterozygosity by Individual",
       x = "Individual",
       y = "Homozygote Observations") +
  scale_y_continuous(labels = comma) 


#OBS_CT (Observation Count): The total number of non-missing loci (i.e. the total number of genotypes evaluated) in the individual.
ggplot(het_data, aes(x = IID_short, y = `OBS_CT`)) + 
  geom_bar(stat = "identity", fill = "plum4") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Total number of non-missing loci per individuals",
       x = "Individual",
       y = "Non-missing loci") +
  scale_y_continuous(labels = comma) 


#F: Inbreeding coefficient for the individual, calculated from homozygous and heterozygous observations. 
#F value suggests an excess of homozygosity (which may indicate inbreeding), 
#while a negative value indicates an excess of heterozygosity over what is expected.
ggplot(het_data, aes(x = IID_short, y = `F`)) + 
  geom_bar(stat = "identity", fill = "#8B8B00") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Inbreeding coefficient per individuals",
       x = "Individual",
       y = "Inbreeding coefficient") +
  scale_y_continuous(labels = comma) 



################################ Proportion of missing genotype in samples after filtering ########################################
###################################################################################################################################

####MISSING_CT: Missing Count. This number indicates how many markers (e.g. SNPs) have missing data for this particular individual. 
#Data is considered missing if the genotype for a given marker could not be determined for the individual.
####OBS_CT: Observation Count. This represents the total number of markers tested or observed for the individual, 
#including those with valid (non-missing) data.
####F_MISS: Fraction Missing. This is the ratio between the number of missing data (MISSING_CT) 
#and the total number of observations (OBS_CT). This value provides a measure of data quality for each individual, 
#with a higher value indicating a greater proportion of missing data.

missing_data <- read.table('missing_filtered_stats.smiss', header = TRUE, sep = "\t")

# Truncate individual IDs to the first 7 characters
missing_data$IID_short <- substr(het_data$IID, 1, 7)

# Create a bar chart of F_MISS for each individual
ggplot(missing_data, aes(x = IID_short, y = F_MISS)) +
  geom_bar(stat = "identity", fill = "steelblue") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) + 
  labs(title = "Fraction of Missing Data by Individual",
       x = "Individual ID",
       y = "Fraction Missing (F_MISS)") + 
  scale_y_continuous(labels = scales::percent_format()) 



