
################# Pre-processing of SweeD result files #####################
######################################################################################

library(dplyr)
library(data.table)

file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/07-outlier_analysis/07.2-SweeD/sub_200k_SNP/SweeD_Report.sub_200k_SNP_west.txt"

# Read the file using readLines to easily manipulate lines
lines <- readLines(file_path)

# Extract scaffolds number and create a list
scaffold_indices <- grep("^//", lines)
scaffold_numbers <- gsub("//", "", lines[scaffold_indices])
data_list <- list()

# Read data without header of before each scaffolds 
for (i in 1:(length(scaffold_indices) - 1)) {
  data_start <- scaffold_indices[i] + 2  
  data_end <- scaffold_indices[i + 1] - 1  
  if (data_start <= data_end) {  
    data_chunk <- read.table(text = lines[data_start:data_end],
                             header = FALSE, comment.char = "/", 
                             col.names = c("Position", "Likelihood", "Alpha", "StartPos", "EndPos"))
    if (nrow(data_chunk) > 0) {  
      data_chunk$scaffold <- scaffold_numbers[i]
      data_list[[i]] <- data_chunk
    } else {
      message(sprintf("The data chunk is empty for the scaffold %s (indices %d à %d).", scaffold_numbers[i], data_start, data_end))
    }
  } else {  
    message(sprintf("No data available for scaffold %s (indices %d à %d).", scaffold_numbers[i], data_start, data_end))
  }
}

last_data_start <- scaffold_indices[length(scaffold_indices)] + 2
if (last_data_start <= length(lines)) {
  # Use tryCatch to handle potential errors while reading
  tryCatch({
    last_data_chunk <- read.table(text = lines[last_data_start:length(lines)],
                                  header = FALSE, comment.char = "/", 
                                  col.names = c("Position", "Likelihood", "Alpha", "StartPos", "EndPos"),
                                  fill = TRUE)  # Add fill = TRUE to handle rows with less than 5 elements
    if (nrow(last_data_chunk) > 0) {
      last_data_chunk$scaffold <- scaffold_numbers[length(scaffold_indices)]
      data_list[[length(scaffold_indices)]] <- last_data_chunk
    } else {
      message(sprintf("The last chunk of data is empty for the scaffold %s.", scaffold_numbers[length(scaffold_indices)]))
    }
  }, error = function(e) {
    message(sprintf("Error reading last data block for scaffold %s : %s", scaffold_numbers[length(scaffold_indices)], e$message))
  })
}

# Combine all data frame
filtered_data <- bind_rows(data_list)

# Filter to keep only lines where Likelihood is different from 0
filtered_data_west <- filtered_data %>% filter(Likelihood != 0)




############################# Graphical analysis ####################################
####################################################################################

library(ggplot2)
library(dplyr)


############################### Graph by scaffold ###################################

# Create the graph for a specific scaffold, for example scaffold 1
scaffold_data <- filtered_data[filtered_data$scaffold == 1,]

ggplot(scaffold_data, aes(x = Position, y = Likelihood)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Likelihood of Selective Sweeps per Position",
       x = "Position along the scaffold",
       y = "Likelihood of Selective Sweep") +
  theme(plot.title = element_text(hjust = 0.5))



############################## Graph on full genome ###############################

# Make sure 'scaffold' is a numeric variable
filtered_data$scaffold <- as.numeric(as.character(filtered_data$scaffold))
str(filtered_data$scaffold)

# Calculate an offset for each scaffold
max_position <- max(filtered_data$Position)
# Create a new column that will contain the adjusted positions
filtered_data$AdjustedPosition <- filtered_data$Position + as.numeric(as.factor(filtered_data$scaffold)) * (max_position + 10000)  # Ajoute un décalage basé sur le numéro du scaffold

# Graph
ggplot(filtered_data, aes(x = AdjustedPosition, y = Likelihood, color = as.factor(scaffold))) +
  geom_point(alpha=0.6) +  
  theme_minimal() +
  labs(x = "Adjusted Position (each scaffold offset)",
       y = "Likelihood",
       color = "Scaffold") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))


####################### Comparing outlier value  ###################################

# West outliers 
unique_west <- anti_join(filtered_data_west, filtered_data_other, by = "Position")
print(unique_west)

# Other sites outliers
unique_other <- anti_join(filtered_data_other, filtered_data_west, by = "Position")
print(unique_other)

# Commun outliers in the 2 data set
common_sweeps <- inner_join(filtered_data_west, filtered_data_other, by = "Position", 
                            suffix = c(".West", ".Other"))
print(common_sweeps)


# Computing west/other ratio
common_sweeps <- common_sweeps %>%
  mutate(Likelihood_Ratio = Likelihood.West / Likelihood.Other)
print(common_sweeps)


















