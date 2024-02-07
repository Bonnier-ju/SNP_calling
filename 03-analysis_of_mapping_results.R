
####Analysis of mapping results####

# Load packages 
library(tidyr)
library(ggplot2)

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








#############################################################
# Analyse des fichiers BED pour la distribution des lectures#
# Cette partie dépend des analyses spécifiques que vous souhaitez effectuer, par exemple, la densité des lectures par région.
# Pour un exemple simple, on pourrait compter le nombre de lectures dans chaque région à partir des fichiers dans statsDir.

bedDir <- "/03-mapping_samples/03-bed_files/"
statsDir <- "/03-mapping_samples/03-bed_analysis/"

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


