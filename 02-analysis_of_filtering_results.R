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

# Read only one log file
process_log_file <- function(file_path) {
  lines <- read_lines(file_path)

  
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

# Apply function in all files and combine results
log_data <- map_df(log_files, process_log_file)

# Data summary
print(log_data)

# Save data in CSV file
write_csv(log_data, file.path(log_dir, "summary_log_data.csv"))

# Plots for passed reads
ggplot(log_data, aes(x = FileName, y = TotalReads_passed, fill = FileName)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() + # Utiliser une échelle de couleurs
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90")) +
  labs(title = "Reads passed by sample",
       x = "Sample",
       y = "Passed reads") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))


#Plot for duplication rates 
ggplot(log_data, aes(x = FileName, y = Duplication_rates, fill = FileName)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
        legend.position = "none", 
        panel.grid.major.y = element_line(color = "black")) +
  labs(title = "Duplication rates (percent)",
       x = "Sample",
       y = "Duplication rates")


#Plot for failed reads
log_data_long <- log_data %>%
  pivot_longer(cols = c(Failedreads_lq, Failedreads_N, Failedread_short),
               names_to = "Measurement",
               values_to = "Value")

ggplot(log_data_long, aes(x = FileName, y = Value, fill = Measurement)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Failed Reads Analysis",
       x = "Sample",
       y = "Count") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))


#Plots for trimmed adapters
log_data_long <- log_data %>%
  pivot_longer(cols = c(Adapters_trimmed_reads, Adapters_trimmed_base),
               names_to = "Measurement",
               values_to = "Value")

ggplot(log_data_long, aes(x = FileName, y = Value, fill = Measurement)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Trimmed adapters",
       x = "Sample",
       y = "Count") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))




