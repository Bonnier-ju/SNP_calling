#####Analysis of filtering results#####

#Installing library 
library(tidyverse)

setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/02-filtering_reads")


# Chemin vers le dossier contenant vos fichiers de log
log_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/02-filtering_reads/log_files_filtered"

# Lister tous les fichiers de log dans le dossier
log_files <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)

# Fonction pour lire et traiter un seul fichier de log
process_log_file <- function(file_path) {
  lines <- read_lines(file_path)
  
  # Exemple de récupération de données, ajustez selon le format de votre fichier log
  total_reads <- str_extract(string = lines[str_detect(lines, "total_reads")], pattern = "\\d+")
  passed_reads <- str_extract(string = lines[str_detect(lines, "passed_filter_reads")], pattern = "\\d+")
  
  data.frame(
    FileName = basename(file_path),
    TotalReads = as.numeric(total_reads),
    PassedReads = as.numeric(passed_reads)
  )
}

# Appliquer la fonction à tous les fichiers et combiner les résultats
log_data <- map_df(log_files, process_log_file)

# Afficher le résumé des données
print(log_data)

# Enregistrer le résumé dans un fichier CSV
write_csv(log_data, file.path(log_dir, "summary_log_data.csv"))

# Analyse et visualisation
# Par exemple, tracer le nombre de lectures passées pour chaque échantillon
ggplot(log_data, aes(x = FileName, y = PassedReads)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Nombre de lectures passées par échantillon",
       x = "Échantillon",
       y = "Lectures passées") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Sauvegarder le graphique
ggsave(file.path(log_dir, "passed_reads_plot.pdf"), width = 10, height = 8)
