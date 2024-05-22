
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
file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.2-Variables_climatiques_guyane/Variable_full_Guyane.csv"
data <- read_csv(file_path, locale = locale(encoding = "latin1"))

# Supprimer les colonnes qui ne sont pas des variables environnementales
data <- data %>%
  select(-wkt_geom, -fid)

head(data)

# Convertir les variables catégorielles en variables numériques (encodage one-hot)
data <- data %>%
  mutate_if(is.character, as.factor) %>%
  mutate_if(is.factor, as.numeric)

# Normaliser les données
data_scaled <- scale(data)

# Réaliser l'ACP
acp_result <- PCA(data, scale.unit = TRUE, ncp = 5, graph = FALSE)
print(acp_result)

# Visualisation des variances expliquées par les composantes principales
fviz_screeplot(acp_result, addlabels = TRUE, ylim = c(0, 50))

# Visualisation des variables sur le plan factoriel 1 et 2
fviz_pca_var(acp_result, axes = c(1, 2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  ggtitle("ACP - Dimensions 1 et 2")

# Visualisation des variables sur le plan factoriel 1 et 3
fviz_pca_var(acp_result, axes = c(1, 3), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  ggtitle("ACP - Dimensions 1 et 3")

# Visualisation des variables sur le plan factoriel 2 et 3
fviz_pca_var(acp_result, axes = c(2, 3), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  ggtitle("ACP - Dimensions 2 et 3")

# Visualisation des variables sur le plan factoriel 2 et 4
fviz_pca_var(acp_result, axes = c(2, 4), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  ggtitle("ACP - Dimensions 2 et 4")

# Visualisation des variables sur le plan factoriel 4 et 5
fviz_pca_var(acp_result, axes = c(4, 5), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  ggtitle("ACP - Dimensions 4 et 5")

# Afficher les contributions des variables
fviz_contrib(acp_result, choice = "var", axes = 1, top = 20)
fviz_contrib(acp_result, choice = "var", axes = 2, top = 20)
fviz_contrib(acp_result, choice = "var", axes = 3, top = 20)
fviz_contrib(acp_result, choice = "var", axes = 4, top = 20)
fviz_contrib(acp_result, choice = "var", axes = 5, top = 20)




