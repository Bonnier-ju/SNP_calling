
####Analysis of mapping results####

# Load packages 
install.packages("ggplot2")
library(ggplot2)

setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/03-mapping_samples")


# Chemin vers les dossiers contenant les fichiers de sortie
flagstatDir <- "/03-mapping_samples/03-flagstat_files/"
bedDir <- "/03-mapping_samples/03-bed_files/"
statsDir <- "/03-mapping_samples/03-bed_analysis/"

# Lire et analyser les fichiers flagstat
# Liste des fichiers flagstat
flagstatFiles <- list.files(flagstatDir, pattern = "*.txt", full.names = TRUE)

# Fonction pour lire un fichier flagstat et extraire des informations clés
readFlagstat <- function(file) {
  lines <- readLines(file)
  totalReads <- as.numeric(strsplit(lines[1], " ")[[1]][1])
  mappedReads <- as.numeric(strsplit(lines[5], " ")[[1]][1])
  data.frame(
    FileName = basename(file),
    TotalReads = totalReads,
    MappedReads = mappedReads,
    MappedPercent = (mappedReads / totalReads) * 100
  )
}

# Appliquer la fonction à tous les fichiers et combiner les résultats
flagstatData <- do.call(rbind, lapply(flagstatFiles, readFlagstat))

# Afficher un résumé
print(flagstatData)

# Visualisation de la proportion de lectures mappées pour chaque échantillon
ggplot(flagstatData, aes(x = FileName, y = MappedPercent)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of Mapped Reads by Sample", x = "Sample", y = "Percentage of Mapped Reads (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Analyse des fichiers BED pour la distribution des lectures
# Cette partie dépend des analyses spécifiques que vous souhaitez effectuer, par exemple, la densité des lectures par région.
# Pour un exemple simple, on pourrait compter le nombre de lectures dans chaque région à partir des fichiers dans statsDir.

# Liste des fichiers de statistiques BED
statsFiles <- list.files(statsDir, pattern = "*_stats.txt", full.names = TRUE)

# Fonction pour lire les statistiques BED
readBedStats <- function(file) {
  data <- read.table(file, header = FALSE)
  colnames(data) <- c("Count", "Region")
  data$FileName <- basename(file)
  return(data)
}

# Appliquer la fonction à tous les fichiers et combiner les résultats
bedStatsData <- do.call(rbind, lapply(statsFiles, readBedStats))

# Afficher un résumé
print(head(bedStatsData))


