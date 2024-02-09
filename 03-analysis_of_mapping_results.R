###################################
####Analysis of mapping results####
###################################


# Load packages 
library(tidyr)
library(ggplot2)
library(dplyr)

setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/03-mapping_samples")


# Chemin vers les dossiers contenant les fichiers de sortie
flagstatDir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/03-mapping_samples/03-flagstat_files/" 

# Lire et analyser les fichiers flagstat
# Liste des fichiers flagstat
flagstatFiles <- list.files(flagstatDir, pattern = "*.txt", full.names = TRUE)
print(flagstatFiles) 


# Fonction pour lire un fichier flagstat et extraire des informations clés
readFlagstat <- function(file) {
  lines <- readLines(file)
  totalReads <- as.numeric(strsplit(lines[1], " ")[[1]][1])
  mappedReads <- as.numeric(strsplit(lines[7], " ")[[1]][1])
  primary_mapped <- as.numeric(strsplit(lines[2], " ")[[1]][1])
  secondary_mapped <- as.numeric(strsplit(lines[3], " ")[[1]][1])
  data.frame(
    FileName = basename(file),
    TotalReads = totalReads,
    MappedReads = mappedReads,
    PrimaryMapped = primary_mapped,
    SecondaryMapped = secondary_mapped,
    MappedPercent = (mappedReads / totalReads) * 100
  )
}

# Appliquer la fonction à tous les fichiers et combiner les résultats
flagstatData <- do.call(rbind, lapply(flagstatFiles, readFlagstat))
flagstatData$FileName <- substr(flagstatData$FileName, 1, 7)
flagstatData$Color <- ifelse(flagstatData$MappedPercent < 80, "Below 80%", "Above 80%")


# Afficher un résumé
print(flagstatData)

###Plot for percentage of mapping reads###
ggplot(flagstatData, aes(x = FileName, y = MappedPercent, fill = Color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Below 80%" = "red", "Above 80%" = "lightblue")) +
  geom_text(aes(label = paste0(sprintf("%.1f", MappedPercent), "%")), 
            position = position_stack(vjust = 0.5), 
            color = "black", 
            angle = 90, 
            size = 4.2) +
  theme_light() +
  labs(title = "Proportion of Mapped Reads by Sample", 
       x = "Sample", 
       y = "Percentage of Mapped Reads (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")) +
  guides(fill = FALSE) 



###Plot for primary and secondary reads###

flagstatData <- do.call(rbind, lapply(flagstatFiles, readFlagstat))
flagstatData$FileName <- substr(flagstatData$FileName, 1, 7)

#Primary mapping 
ggplot(flagstatData, aes(x = FileName, y = PrimaryMapped)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Primary Mapped Reads by Sample", 
       x = "Sample", 
       y = "Number of Primary Mapped Reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Secondary mapping
#alignements supplémentaires pour des lectures qui peuvent s'aligner à plusieurs endroits du génome#

ggplot(flagstatData, aes(x = FileName, y = SecondaryMapped)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_minimal() +
  labs(title = "Secondary Mapped Reads by Sample", 
       x = "Sample", 
       y = "Number of Secondary Mapped Reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






###################################################################
#### Analyse des fichiers BED pour la distribution des lectures####


# Définir le chemin vers le répertoire contenant les fichiers de statistiques
stats_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/03-mapping_samples/03-bed_analysis/"

# Lister tous les fichiers de statistiques dans le répertoire
stats_files <- list.files(stats_dir, pattern = "*.txt", full.names = TRUE)

# Llire et préparer les données d'un fichier de statistiques
read_stats_file <- function(file_path) {
  data <- read.table(file_path, col.names = c("ReadCount", "Region"))
  data$Sample <- gsub(".*/|_stats\\.txt$", "", file_path)  # 
  return(data)
}

# Lire les données de tous les fichiers et les combiner en un seul dataframe
all_data <- do.call(rbind, lapply(stats_files, read_stats_file))

# Analyse des données

## Résumé statistique des comptages de lectures par région
summary_stats <- all_data %>%
  group_by(Region) %>%
  summarise(MeanReadCount = mean(ReadCount), 
            MedianReadCount = median(ReadCount), 
            SDReadCount = sd(ReadCount),
            MinReadCount = min(ReadCount),
            MaxReadCount = max(ReadCount))

print(summary_stats)

## Visualisation de la distribution des comptages de lectures par région
ggplot(all_data, aes(x = Region, y = ReadCount)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution des comptages de lectures par région",
       x = "Région",
       y = "Comptage de lectures") +
  scale_y_log10()  # Utiliser l'échelle logarithmique pour les comptages de lectures


## Comparaison des comptages de lectures entre les échantillons
ggplot(all_data, aes(x = Sample, y = ReadCount, fill = Region)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Comparaison des comptages de lectures entre les échantillons",
       x = "Échantillon",
       y = "Comptage de lectures") +
  facet_wrap(~Region, scales = "free_y")  # Séparer les graphiques par région







