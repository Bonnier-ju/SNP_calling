########################################################################################################
################################# Population structure analysis ########################################
########################################################################################################



################################## PCA analysis for individuals ###############################
###############################################################################################

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)
library(data.table)
library(scales) 
library(ggrepel)
library(tidyverse)

# Define the base path for input files
setwd("C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.2-global_population_stats/06.2.1-PCA")

base_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.2-global_population_stats/06.2.1-PCA"

# Define the path for the population information CSV file
population_info_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/Pop_file.csv"


# Read PCA results and population information
pca_data <- read_delim(file.path(base_path, "sub_pruned_snp.eigenvec"), col_names = FALSE, 
                       delim = " ") # Read PCA results with a tab delimiter
population_info <- read_csv(population_info_path) # Read population info, assuming CSV format

# Prepare PCA data for merging by renaming IID to SampleID
# and creating a new column with the first 7 letters of each SampleID
pca_data <- pca_data %>%
  rename(SampleID = X2) %>%
  mutate(SampleID_short = substr(SampleID, 7, 9)) # Extract characters 7 to 9 from SampleID


# Merge PCA data with population information based on SampleID
pca_data <- pca_data %>%
  inner_join(population_info, by = "SampleID_short")

pca_data <- pca_data %>%
  rename(PC1 = X3, PC2 = X4, PC3 = X5,)


pop_order <- c("Apatou", "Acarouany", "Piste_St_Elie", "Nouragues_Inselberg", "Cacao", "Regina", "Saut_Lavilette", "Foret_Regina_St_Georges", "MC_87", "MC_88", "St_georges")

pca_data$Pop <- factor(pca_data$Site, levels = pop_order)


my_colors <- c("Apatou" = "#8B2500", "Acarouany" = "#CD3700",  "Piste_St_Elie" = "#CDCD00",
               "Nouragues_Inselberg" = "green4", "MC_88" = "#A020F0","Cacao" = "#838B8B",
               "MC_87" = "#DDA0DD", "Foret_Regina_St_Georges" = "deepskyblue3","Regina" = "#0000FF",
               "Saut_Lavilette" = "#96CDCD",  "St_georges" = "#551A8B")


#Variance %
eigenvalues <- read.table("sub_pruned_snp.eigenval", header = FALSE)
eigenvalues <- eigenvalues$V1

total_variance <- sum(eigenvalues)
percent_variance <- (eigenvalues / total_variance) * 100

pc1_percent <- percent_variance[1]
pc2_percent <- percent_variance[2]
pc3_percent <- percent_variance[3]


# Visualisation PCA avec ggplot2, en utilisant la palette de couleurs personnalisée
ggplot(pca_data, aes(x = PC2, y = PC3, label = SampleID_short, color = Pop)) +
  geom_point() + 
  geom_text_repel(size = 3) +
  #stat_ellipse(type = "t", level = 0.95) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "PCA Plot", 
       x = paste("PC2 (", round(pc2_percent, 2), "%)", sep=""), 
       y = paste("PC3 (", round(pc3_percent, 2), "%)", sep=""),
       color = "Population") + 
  scale_color_manual(values = my_colors) + # Utiliser les couleurs définies manuellement
  theme(legend.position = "right", legend.title = element_text(face = "bold"))



############################## Find the optimal number of clusters #####################################
########################################################################################################
#Script following Fishery genomics of Sebastes: investigating population structure using ADMIXTURE by laurabenestan
# https://github.com/laurabenestan/Admixture


library(stringr)
library(ggplot2)
library(dplyr)

cv <- read.table("cross_validation.txt")


#Analyze the cross-validation results Then, add a K-cluster column indicating the number of K you test 
#and select only two columns of interest, CV and K.
cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
CV <- select(cv, V4,K)

#Rename your two columns CV and K-cluster
colnames(CV) <- c("CV","K")

#Do a graph showing the cross validation results. Then select the optimal number of clusters regarding :
# - the lowest cross validation error
# - when the cross-validation error decrease the most


# Convert 'K' to numeric, removing any non-numeric prefix if necessary
CV$K <- as.numeric(gsub("K=", "", CV$K))

# Définir les titres
graph_title <- "Cross-validation Error vs. Number of Clusters"
x_title <- "Number of Clusters (K)"
y_title <- "Cross-validation Error"

# Recreate the plot
graph_1 <- ggplot(CV, aes(x = K, y = CV)) +
  geom_line(size = 1) +  # Specify line width
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +  # Ensure 'K' is now numeric
  labs(title = graph_title, x = x_title, y = y_title) +
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        legend.title = element_blank())

# Print the plot
print(graph_1)

############################## Clustering visualization ################################################
########################################################################################################


library(ggplot2)
library(dplyr)
library(tidyr)

# Chemins vers les résultats d'admixture, les informations sur la population, et le fichier des sites
admixture_results_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.1-admixture/sub_pruned_500000SNP"
pop_files_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/Pop_file.csv"

# Lecture du fichier des sites d'échantillonnage
sites_data <- read.csv(pop_files_path, header = TRUE, stringsAsFactors = FALSE)
colnames(sites_data) <- c("SamplingSite", "IndividualID")

# Choisir une valeur spécifique de K pour la visualisation
K <- 3

# Lecture des résultats d'admixture avec les noms des individus
admixture_file <- sprintf("%s.%d.named.Q", admixture_results_path, K)
admixture_data <- read.table(admixture_file, col.names = c("IndividualID", paste0("Ancestry", 1:K)))

# Modification des identifiants des individus pour ne garder que les caractères 7, 8 et 9
admixture_data$IndividualID <- substr(admixture_data$IndividualID, 7, 9)

# Fusion des données d'admixture avec les sites d'échantillonnage
admixture_data <- merge(admixture_data, sites_data, by = "IndividualID")

# Transformation des données pour la visualisation
admixture_long <- admixture_data %>%
  gather(key = "Ancestry", value = "Proportion", -IndividualID, -SamplingSite) %>%
  mutate(Ancestry = factor(Ancestry, levels = paste0("Ancestry", 1:K)))

# Filtre pour exclure les données où SamplingSite est "AFD"
admixture_long <- admixture_long %>% 
  filter(IndividualID != "AFD")

# Assurez-vous que vos données sont dans admixture_long et que pop_order est défini
pop_order <- c("Apatou", "Acarouany",
               "Saut_Lavilette", "Foret_Regina_St_Georges", "MC_87", "MC_88", "St_georges", "Piste_St_Elie", "Nouragues_Inselberg", "Cacao", "Regina")

# Réordonner les niveaux du facteur SamplingSite selon pop_order et ordonner les données
admixture_long$SamplingSite <- factor(admixture_long$SamplingSite, levels = pop_order)
admixture_long <- admixture_long %>% arrange(SamplingSite, IndividualID)

# Création du graphique en utilisant ggplot
ggplot(admixture_long, aes(x = reorder(IndividualID, as.numeric(SamplingSite)), y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),  # Rotation des étiquettes pour une meilleure lisibilité
        axis.title.x = element_blank(),  # Retirer le titre de l'axe des x
        axis.ticks.x = element_blank(),  # Retirer les ticks de l'axe des x
        panel.grid.major.x = element_blank(),  # Retirer les grilles majeures en x
        panel.grid.minor.x = element_blank(),  # Retirer les grilles mineures en x
        legend.position = "bottom") +  # Position de la légende
  labs(y = "Proportion d'ancestralité", fill = "Ancestralité")  # Étiquettes de l'axe des y et de la légende



################################### Plot admixture results on map #######################################
#########################################################################################################

# => échec 


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggmap)
library(ggrepel)

pop_files_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/geoloc_site.csv"

geoloc_site <- read.csv(pop_files_path, stringsAsFactors = FALSE)
colnames(geoloc_site) <- c('Site', 'Latitude', 'Longitude')

# Calculer les moyennes d'ancestralité pour chaque site et chaque composante d'ancestralité
ancestry_averages <- admixture_long %>%
  group_by(SamplingSite, Ancestry) %>%
  summarise(MeanProportion = mean(Proportion, na.rm = TRUE)) %>%
  ungroup()

# Fusionner les moyennes avec les coordonnées GPS
ancestry_averages <- merge(ancestry_averages, geoloc_site, by.x = "SamplingSite", by.y = "Site")

register_google(key = "AIzaSyA2DS11azaPzP0gT2UjRCNi42-H_xcGA50")

# Obtenir une carte de la Guyane Française
map <- get_map(location = 'French Guiana', zoom = 7, maptype = "terrain")

ggmap(map) +
  geom_point(data = ancestry_averages, 
             aes(x = Longitude, y = Latitude), 
             color = "blue", size = 2, alpha = 0.7) +  # Points de taille fixe
  theme_minimal() +
  theme(legend.position = "none")  # Aucune légende pour les points
































