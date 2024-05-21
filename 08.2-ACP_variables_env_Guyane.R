
################# Réalisation d'une ACP sur les variables climatiques en Guyane #####################
#####################################################################################################

install.packages("FactoMineR")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("readr")
install.packages("dplyr")
install.packages("ggplot2")

library(readr)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)

# Lire le fichier CSV
file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/Variables_climatiques_guyane/Variable_full_Guyane.csv"
data <- read_csv(file_path, locale = locale(encoding = "latin1"))

# Supprimer les colonnes qui ne sont pas des variables environnementales
data <- data %>%
  select(-wkt_geom, -fid)

# Vérifier les premières lignes du dataframe pour s'assurer que les données ont été lues correctement
head(data)

# Réaliser l'ACP
acp_result <- PCA(data, scale.unit = TRUE, ncp = 5, graph = FALSE)

# Résumer les résultats de l'ACP
print(acp_result)

# Visualisation des variances expliquées par les composantes principales
fviz_screeplot(acp_result, addlabels = TRUE, ylim = c(0, 50))

# Visualisation des variables sur le plan factoriel
fviz_pca_var(acp_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# Visualisation des individus sur le plan factoriel
fviz_pca_ind(acp_result, geom.ind = "point", col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# Afficher les contributions des variables
fviz_contrib(acp_result, choice = "var", axes = 1, top = 20)
fviz_contrib(acp_result, choice = "var", axes = 2, top = 20)



#########################################################################################
# Charger les bibliothèques nécessaires
library(readr)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)

# Lire les données à partir du fichier CSV
file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.2-Variables_climatiques_guyane/Variable_full_Guyane.csv"
data <- read_csv(file_path)

# Afficher les premières lignes des données pour comprendre leur structure
print(head(data))

# Convertir les variables catégorielles en variables numériques (encodage one-hot)
data <- data %>%
  mutate_if(is.character, as.factor) %>%
  mutate_if(is.factor, as.numeric)

# Normaliser les données
data_scaled <- scale(data)

# Réaliser l'ACP
pca_result <- PCA(data_scaled, graph = FALSE)

# Visualiser les résultats de l'ACP
# Biplot des individus et des variables
fviz_pca_biplot(pca_result, repel = TRUE, col.var = "blue", col.ind = "red") +
  ggtitle("Biplot des individus et des variables")

# Variance expliquée par chaque composante principale
fviz_screeplot(pca_result, addlabels = TRUE, ylim = c(0, 50)) +
  ggtitle("Variance expliquée par chaque composante principale")

# Contributions des variables à la première composante principale
fviz_pca_var(pca_result, col.var = "contrib") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red", midpoint = 50) +
  theme_minimal() +
  ggtitle("Contributions des variables à la première composante principale")

# Contributions des individus à la première composante principale
fviz_pca_ind(pca_result, col.ind = "cos2") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red", midpoint = 0.5) +
  theme_minimal() +
  ggtitle("Contributions des individus à la première composante principale")



