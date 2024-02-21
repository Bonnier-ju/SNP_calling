#####Analysis of filtering results#####

#Installing library 
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)

setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/02-filtering_reads")


# Path of log files from fastq results 
log_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/02-filtering_reads/log_files_filtered"

# All log files in the folder
log_files <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)


process_log_file <- function(file_path) {
  lines <- readr::read_lines(file_path)
  
  # Extraction des données
  total_reads_passed <- str_extract(string = lines[str_detect(lines, "reads passed filter")], pattern = "\\d+")
  failed_read_lq <- str_extract(string = lines[str_detect(lines, "reads failed due to low quality")], pattern = "\\d+")
  failed_read_N <- str_extract(string = lines[str_detect(lines, "reads failed due to too many N")], pattern = "\\d+")
  failed_read_short <- str_extract(string = lines[str_detect(lines, "reads failed due to too short")], pattern = "\\d+")
  reads_adapters_trimmed <- str_extract(string = lines[str_detect(lines, "reads with adapter trimmed")], pattern = "\\d+")
  based_adapters_trimmed <- str_extract(string = lines[str_detect(lines, "bases trimmed due to adapters")], pattern = "\\d+")
  duplication_rates <- str_extract(string = lines[str_detect(lines, "Duplication rate")], pattern = "\\d+\\.?\\d*")
  
  # Vérification et correction pour s'assurer que toutes les variables ont une longueur de 1
  variables <- list(total_reads_passed, failed_read_lq, failed_read_N, failed_read_short, reads_adapters_trimmed, based_adapters_trimmed, duplication_rates)
  variables <- lapply(variables, function(x) if(is.na(x) || length(x) == 0) NA else as.numeric(x))
  
  data.frame(
    FileName = basename(file_path),
    TotalReads_passed = variables[[1]],
    Failedreads_lq = variables[[2]],
    Failedreads_N = variables[[3]],
    Failedread_short = variables[[4]],
    Adapters_trimmed_reads = variables[[5]],
    Adapters_trimmed_base = variables[[6]],
    Duplication_rates = variables[[7]]
  )
}


# Apply function in all files and combine results
log_data <- map_df(log_files, process_log_file)

log_data <- log_data %>%
  mutate(ShortFileName = substr(FileName, 1, 7))


# Data summary
print(log_data)

# Save data in CSV file
write_csv(log_data, file.path(log_dir, "summary_log_data.csv"))

# Plots for passed reads
#Read Passed Filter: This term refers to sequence reads that have passed the quality and filtering criteria specific to the sequencing platform. 
#Reads that "pass the filter" are considered to be of good quality and suitable for further analysis.

ggplot(log_data, aes(x = ShortFileName, y = TotalReads_passed, fill = FileName)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() + # Utiliser une échelle de couleurs
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 11),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90")) +
  labs(title = "Reads passed by sample",
       x = "Sample",
       y = "Passed reads") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))


#Plot for duplication rates 
ggplot(log_data, aes(x = ShortFileName, y = Duplication_rates, fill = FileName)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 11),
        legend.position = "none", 
        panel.grid.major.y = element_line(color = "black")) +
  labs(title = "Duplication rates (percent)",
       x = "Sample",
       y = "Duplication rates")


#Plot for failed reads
#Read Failed Due to Low Quality: These reads have failed to meet the minimum quality criteria, often due to an average quality score that is too low over the length of the read. 
#Quality scores are generally assigned to each base and indicate confidence in the accuracy of that base.

#Reads Failed Due to Too Many N: Reads containing an excessive number of "N" bases fail filtering. 
#An "N" base indicates uncertainty as to the identity of the base (A, T, C, or G), often due to a weak or ambiguous signal during sequencing.

#Reads Failed Due to Too Short: These reads are rejected because they are shorter than the minimum length defined for the analysis. 
#The minimum length is often specified to ensure that reads are sufficiently informative for reliable alignment or assembly.
log_data_long <- log_data %>%
  pivot_longer(cols = c(Failedreads_lq, Failedreads_N, Failedread_short),
               names_to = "Measurement",
               values_to = "Value")

ggplot(log_data_long, aes(x = ShortFileName, y = Value, fill = Measurement)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Failed Reads Analysis",
       x = "Sample",
       y = "Count") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))


#Plots for trimmed adapters
ggplot(log_data, aes(x = ShortFileName, y = Adapters_trimmed_reads, fill = FileName)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() + # Utiliser une échelle de couleurs
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 11),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90")) +
  labs(title = "Adaptaters trimmed",
       x = "Sample",
       y = "Nb Adapters") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))




