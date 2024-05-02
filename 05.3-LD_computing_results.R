
# Charger les bibliothèques nécessaires
library(ggplot2)

file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.4-LD_pruning/ld_results_10percent.txt"

# Lire les données depuis un fichier CSV ou TXT

data <- read.table(file_path, header = FALSE, sep = "", skip = 1,
                   col.names = c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"), 
                   fill = TRUE, strip.white = TRUE)

# Afficher les premières lignes pour s'assurer que les données sont bien chargées
head(data)

# Calculer la distance entre les SNP
data$Distance <- abs(data$BP_A - data$BP_B)

# Filtrer pour conserver uniquement les distances inférieures ou égales à 10000
data <- data[data$Distance <= 10000, ]

summary(data$Distance)

# Sous-échantillonnage d'un quart des données
sample_size <- nrow(data) * 0.25
sample_indices <- sample(1:nrow(data), size = sample_size)
sampled_data <- data[sample_indices, ]

# Créer un graphique de dispersion de R^2 en fonction de la distance
ggplot(sampled_data, aes(x = Distance, y = R2)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +  
  labs(x = "Distance entre SNP (bp)", y = "Valeur de R²",
       title = "Relation entre le déséquilibre de liaison (R²) et la distance entre SNP") +
  theme(plot.title = element_text(hjust = 0.5))
