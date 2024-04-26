
################# Pré-traitement des fichiers resultats de SweeD #####################
######################################################################################

# Charger les bibliothèques nécessaires
library(dplyr)
library(data.table)

# Chemin vers le fichier de résultats de SweeD
file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/07-outlier_analysis/07.2-SweeD/sub_200k_SNP/SweeD_Report.sub_200k_SNP_west.txt"

# Lire le fichier en utilisant readLines pour manipuler facilement les lignes
lines <- readLines(file_path)

# Identifier les indices des lignes contenant les identifiants de scaffold
scaffold_indices <- grep("^//", lines)

# Extraire juste les numéros des scaffolds
scaffold_numbers <- gsub("//", "", lines[scaffold_indices])

# Créer une liste pour stocker les données par scaffold
data_list <- list()

# Lire les données entre chaque index de scaffold, en utilisant read.table avec comment.char = "/"
for (i in 1:(length(scaffold_indices) - 1)) {
  data_start <- scaffold_indices[i] + 2  # Début des données après la ligne de noms de colonnes
  data_end <- scaffold_indices[i + 1] - 1  # Fin des données juste avant le prochain identifiant de scaffold
  if (data_start <= data_end) {  # Vérifier qu'il y a des lignes de données à lire
    data_chunk <- read.table(text = lines[data_start:data_end],
                             header = FALSE, comment.char = "/", 
                             col.names = c("Position", "Likelihood", "Alpha", "StartPos", "EndPos"))
    if (nrow(data_chunk) > 0) {  # S'assurer que le data_chunk contient des données
      data_chunk$scaffold <- scaffold_numbers[i]
      data_list[[i]] <- data_chunk
    } else {
      message(sprintf("Le chunk de données est vide pour le scaffold %s (indices %d à %d).", scaffold_numbers[i], data_start, data_end))
    }
  } else {  # Gérer le cas où aucun data n'est disponible entre les indices
    message(sprintf("Aucune donnée disponible pour le scaffold %s (indices %d à %d).", scaffold_numbers[i], data_start, data_end))
  }
}

# Ajouter le dernier chunk après le dernier index de scaffold
last_data_start <- scaffold_indices[length(scaffold_indices)] + 2
if (last_data_start <= length(lines)) {
  # Utiliser tryCatch pour gérer les erreurs potentielles lors de la lecture
  tryCatch({
    last_data_chunk <- read.table(text = lines[last_data_start:length(lines)],
                                  header = FALSE, comment.char = "/", 
                                  col.names = c("Position", "Likelihood", "Alpha", "StartPos", "EndPos"),
                                  fill = TRUE)  # Ajouter fill = TRUE pour gérer les lignes avec moins de 5 éléments
    if (nrow(last_data_chunk) > 0) {
      last_data_chunk$scaffold <- scaffold_numbers[length(scaffold_indices)]
      data_list[[length(scaffold_indices)]] <- last_data_chunk
    } else {
      message(sprintf("Le dernier chunk de données est vide pour le scaffold %s.", scaffold_numbers[length(scaffold_indices)]))
    }
  }, error = function(e) {
    message(sprintf("Erreur lors de la lecture du dernier bloc de données pour le scaffold %s : %s", scaffold_numbers[length(scaffold_indices)], e$message))
  })
}

# Combiner tous les dataframes dans un seul dataframe
filtered_data <- bind_rows(data_list)

# Filtrer pour garder seulement les lignes où Likelihood est différent de 0
filtered_data_west <- filtered_data %>% filter(Likelihood != 0)

# Vérifier et sauvegarder les résultats
print(head(filtered_data))
write.csv(filtered_data, "FilteredDataByScaffold.csv", row.names = FALSE)




############################# Analyse Graphique ####################################
####################################################################################

library(ggplot2)
library(dplyr)


#### Graph by scaffold

# Créer le graphique pour un scaffold spécifique, par exemple scaffold 1
scaffold_data <- filtered_data[filtered_data$scaffold == 1,]

# Générer le graphique
p <- ggplot(scaffold_data, aes(x = Position, y = Likelihood)) +
  geom_point() +  # Ajouter des points pour visualiser les données
  theme_minimal() +  # Utiliser un thème minimal pour une meilleure visibilité
  labs(title = "Likelihood of Selective Sweeps per Position",
       x = "Position along the scaffold",
       y = "Likelihood of Selective Sweep") +
  theme(plot.title = element_text(hjust = 0.5))  # Centrer le titre du graphique

# Afficher le graphique
print(p)




### Graph of full Genome

# Assurez-vous que 'scaffold' est une variable numérique
filtered_data$scaffold <- as.numeric(as.character(filtered_data$scaffold))

# Vérifiez que la conversion est réussie
str(filtered_data$scaffold)

# Calculer un décalage pour chaque scaffold
max_position <- max(filtered_data$Position)
# Créer une nouvelle colonne qui contiendra les positions ajustées
filtered_data$AdjustedPosition <- filtered_data$Position + as.numeric(as.factor(filtered_data$scaffold)) * (max_position + 10000)  # Ajoute un décalage basé sur le numéro du scaffold

# Générer le graphique avec des positions décalées
ggplot(filtered_data, aes(x = AdjustedPosition, y = Likelihood, color = as.factor(scaffold))) +
  geom_point(alpha=0.6) +  
  theme_minimal() +
  labs(x = "Adjusted Position (each scaffold offset)",
       y = "Likelihood",
       color = "Scaffold") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))


########## Comparing outlier value in west and other ###########

# Sweeps uniques à l'ouest
unique_west <- anti_join(filtered_data_west, filtered_data_other, by = "Position")

# Afficher les résultats
print(unique_west)


# Sweeps uniques aux autres régions
unique_other <- anti_join(filtered_data_other, filtered_data_west, by = "Position")

# Afficher les résultats
print(unique_other)

# Sweeps communs
common_sweeps <- inner_join(filtered_data_west, filtered_data_other, by = "Position", 
                            suffix = c(".West", ".Other"))
# Afficher les résultats
print(common_sweeps)


# Ajouter une nouvelle colonne pour le rapport des likelihoods
common_sweeps <- common_sweeps %>%
  mutate(Likelihood_Ratio = Likelihood.West / Likelihood.Other)

# Afficher le résultat
print(common_sweeps)


















