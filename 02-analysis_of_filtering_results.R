#####Analysis of filtering results#####

#Installing library 
library(tidyverse)

setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/02-filtering_reads")


# Path of log files from fastq results 
log_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/02-filtering_reads/log_files_filtered"

# All log files in the folder
log_files <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)

# Read only one log file
process_log_file <- function(file_path) {
  lines <- read_lines(file_path)
  
  
  # S'assurer que `lines` est bien un vecteur de chaînes de caractères
  if (!is.character(lines)) {
    stop("Expected 'lines' to be a character vector, but it was not.")
  }
  

  
  # Data extraction
  total_reads_passed <- str_extract(string = lines[str_detect(lines, "reads passed filter")], pattern = "\\d+")
  failed_read_lq <- str_extract(string = lines[str_detect(lines, "reads failed due to low quality")], pattern = "\\d+")
  failed_read_N <- str_extract(string = lines[str_detect(lines, "reads failed due to too many N")], pattern = "\\d+")
  failed_read_short <- str_extract(string = lines[str_detect(lines, "reads failed due to too short")], pattern = "\\d+")
  reads_adapters_trimmed <- str_extract(string = lines[str_detect(lines, "reads with adapter trimmed")], pattern = "\\d+")
  based_adapters_trimmed <- str_extract(string = lines[str_detect(lines, "bases trimmed due to adapters")], pattern = "\\d+")
  duplication_rates <- str_extract(string = lines[str_detect(lines, "Duplication rate")], pattern = "\\d+")
 
  
  data.frame(
    FileName = basename(file_path),
    TotalReads_passed = as.numeric(total_reads_passed),
    Failedreads_lq = as.numeric(failed_read_lq),
    Failedreads_N = as.numeric(failed_read_N),
    Failedread_short = as.numeric(failed_read_short),
    Adapters_trimmed_reads = as.numeric(reads_adapters_trimmed),
    Adapters_trimmed_base = as.numeric(based_adapters_trimmed),
    Duplication_rates = as.numeric(duplication_rates)
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
ggplot(log_data, aes(x = FileName, y =TotalReads_passed)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Nombre de lectures passées par échantillon",
       x = "Échantillon",
       y = "Lectures passées") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


