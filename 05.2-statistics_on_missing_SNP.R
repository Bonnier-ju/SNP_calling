# Define the base path for input files
base_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.2-normal_filtering"

# Define the path for the population information CSV file
population_info_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.2-normal_filtering/Population_info.csv"


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)
library(data.table)


# PCA analysis for individuals

# Lecture du fichier .eigenvecet csv        
pca_data <- read_delim(file.path(base_path, "final_filtered_filtered_stats.eigenvec"), 
                       delim = "\t")         
         
population_info <- read_csv(population_info_path)

# Préparation des données PCA pour la jointure
pca_data <- pca_data %>%
  rename(SampleID = IID) # Assurez-vous que cette colonne correspond exactement à celle dans `population_info`

# Jointure des données PCA avec les informations de population
pca_data <- pca_data %>%
  inner_join(population_info, by = "SampleID")


# Visualisation PCA avec ggplot2
ggplot(pca_data, aes(x = PC3, y = PC1, color = Pop)) +
  geom_point() + 
  stat_ellipse(type = "t", level = 0.95) + # Ajoute une ellipse de confiance autour des groupes
  geom_hline(yintercept = 0, linetype = "dashed") + # Ajoute une ligne horizontale à y = 0
  geom_vline(xintercept = 0, linetype = "dashed") + # Ajoute une ligne verticale à x = 0
  theme_minimal() +
  labs(title = "PCA Plot", 
       x = "PC1", 
       y = "PC2", 
       color = "Population") + # Ceci assure que la légende est intitulée "Population"
  scale_color_viridis_d(begin = 0.1, end = 0.9, name = "Population", guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.position = "right", legend.title = element_text(face = "bold")) # Personnalise la légende

# Pour sauvegarder le graphique
ggsave("PCA_plot.png", width = 10, height = 8, dpi = 300)

