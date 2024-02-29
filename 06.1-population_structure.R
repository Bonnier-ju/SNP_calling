########################################################################################################
################################# Population structure analysis ########################################
########################################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/")

samples_sites <- read.csv('Population_info.csv', header = TRUE, stringsAsFactors = FALSE)

# Choisir une valeur spécifique de K pour la visualisation
K <- 3 # Ajustez cette valeur en fonction de l'analyse à visualiser

admixture_file <- sprintf("%s.%d.Q", admixture_results_path, K)

# Read admixture results 
admixture_data <- read.table(admixture_file, col.names = paste0("Ancestry", 1:K))

admixture_data <- cbind(IndividualID = samples_sites$ID, admixture_data)

# Link data and sampling sites 
admixture_data <- merge(admixture_data, samples_sites, by.x = "IndividualID", by.y = "ID")

# Translate data to visualisation 
admixture_long <- admixture_data %>%
  gather(key = "Ancestry", value = "Proportion", -IndividualID, -Site) %>%
  mutate(Ancestry = factor(Ancestry, levels = paste0("Ancestry", 1:K)))

# Plot
ggplot(admixture_long, aes(x = IndividualID, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Individu", y = "Proportion d'ancestralité", fill = "Ancestralité")


