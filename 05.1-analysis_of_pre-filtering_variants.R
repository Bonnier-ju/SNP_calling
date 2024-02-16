##################################################################
################# Analysis of pre-filtering results ##################
##################################################################


library(tidyverse)
library(ggplot2)


################# Pre_filtering of variants ######################

# Définir le chemin vers le dossier contenant les fichiers de sortie des taitements sur le clusters
output_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.1-pre_filtering/Input_files/"

# Charger les données
missing_indv <- read_tsv(paste0(output_dir, "missing_indv.imiss.txt"), col_names = TRUE)
missing_site <- read_tsv(paste0(output_dir, "missing_site.lmiss"), col_names = TRUE)
site_quality <- read_tsv(paste0(output_dir, "site_quality.lqual"), col_names = TRUE)
info_fields <- read_tsv(paste0(output_dir, "INFO_fields.INFO.txt"), col_names = TRUE)

# Visualisation du taux de données manquantes par individu
ggplot(missing_indv, aes(x = N_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Taux de données manquantes par individu",
       x = "Taux de données manquantes",
       y = "Nombre d'individus") +
  theme_minimal()

# Visualisation du taux de données manquantes par site
ggplot(missing_site, aes(x = N_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "pink", color = "black") +
  labs(title = "Taux de données manquantes par site",
       x = "Taux de données manquantes",
       y = "Nombre de sites") +
  theme_minimal()

# Visualisation de la qualité des sites
ggplot(site_quality, aes(x = QUAL)) +
  geom_histogram(binwidth = 1, fill = "lightgreen", color = "black") +
  labs(title = "Qualité des sites",
       x = "Qualité",
       y = "Nombre de sites") +
  theme_minimal()

##### Visualisation des champs INFO #####

info_fields_long <- info_fields %>%
  gather(key = "Metric", value = "Value", AC, AF, QD, FS, SOR)

info_plot <- info_fields_long %>%
  ggplot(aes(x = Value, fill = Metric)) +
  geom_histogram(binwidth = 0.1, color = "black") +
  scale_fill_manual(values = c(AC = "#98FB98", AF = "aquamarine3", QD = "#458B74", FS = "#9BCD9B", SOR = "olivedrab3")) +
  facet_wrap(~ Metric, scales = "free_x", labeller = labeller(Metric = c(AC = "Allele Count",
                                                                         AF = "Allele Frequency",
                                                                         QD = "Quality by Depth",
                                                                         FS = "Fisher Strand",
                                                                         SOR = "Strand Odds Ratio"))) +
  labs(title = "Distribution des métriques",
       x = "Valeur",
       y = "Fréquence") +
  theme_minimal() +
  theme(legend.position = "none") # Cache la légende

print(info_plot)

#AC (Allele Count) : Le nombre d'occurrences d'un allèle alternatif dans les échantillons analysés. 
#Pour les variants bialléliques, cela représente le nombre total d'occurrences de l'allèle non référence.
#
#AF (Allele Frequency) : La fréquence de l'allèle alternatif dans l'ensemble des échantillons. 
#Cela donne une idée de la prévalence de l'allèle variant par rapport à l'allèle de référence dans la population étudiée.
#
#QD (Quality by Depth) : La qualité du variant divisée par la profondeur de couverture. 
#Estimation de la confiance dans l'appel du variant en tenant compte de la profondeur de séquençage. 
#Des valeurs élevées indiquent généralement une plus grande confiance.
#
#FS (Fisher Strand) : Le biais de strand, évalué à l'aide d'un test exact de Fisher. 
#Mesure si les allèles variants sont distribués de manière égale entre les brins forward et reverse du séquençage. 
#Des valeurs élevées peuvent indiquer un biais de séquençage, ce qui peut réduire la confiance dans l'appel du variant.
#
#SOR (Strand Odds Ratio) : Un autre indicateur du biais de strand, calculé comme le rapport des chances (odds ratio) 
#de biais entre les brins forward et reverse. Comme pour le FS, des valeurs élevées peuvent suggérer 
#un biais dans la distribution des allèles entre les brins, ce qui peut affecter la fiabilité de l'appel du variant.