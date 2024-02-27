##################################################################
################# Analysis of pre-filtering results ##################
##################################################################


library(tidyverse)
library(ggplot2)
library(dplyr)


################# Pre_filtering of variants ######################

# Définir le chemin vers le dossier contenant les fichiers de sortie des taitements sur le clusters
output_dir <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.1-pre_filtering/Input_files/"

# Charger les données
missing_indv <- read_tsv(paste0(output_dir, "missing_indv.imiss"), col_names = TRUE)
missing_site <- read_tsv(paste0(output_dir, "missing_site.lmiss"), col_names = TRUE)
site_quality <- read_tsv(paste0(output_dir, "site_quality.lqual"), col_names = TRUE)
info_fields <- read_tsv(paste0(output_dir, "INFO_fields.INFO"), col_names = TRUE)

# Visualisation du taux de données manquantes par individu
ggplot(missing_indv, aes(x = N_MISS)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(title = "Missing data per individuals",
       x = "Missing data",
       y = "Number of individuals") +
  theme_minimal()


# Calculate the average quality by sample
quality_avg_by_scaffold <- site_quality %>%
  group_by(CHROM) %>%
  summarise(Mean_QUAL = mean(QUAL, na.rm = TRUE)) # Calculate the average, excluding missing values

# Visualize the average quality by sample
ggplot(quality_avg_by_scaffold, aes(x = CHROM, y = Mean_QUAL)) +
  geom_bar(stat = "identity", fill = "cornflowerblue", color = "black") +
  labs(title = "Average Quality by Scaffold",
       x = "Scaffold",
       y = "Average Quality") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Angle x-axis labels for better readability



##### Visualisation des champs INFO #####

info_fields_long <- info_fields %>%
  gather(key = "Metric", value = "Value", AC, AF, QD, FS, SOR)

info_plot <- info_fields_long %>%
  ggplot(aes(x = Value, fill = Metric)) +
  geom_histogram(binwidth = 30, color = "black") +
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