
#####Analysis of variants calling results#####

# Installation du package si nécessaire
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
BiocManager::install("ggbio")
install.packages("ggplot2")
install.packages("dplyr")

library(VariantAnnotation)
library(ggbio)
library(ggplot2)
library(dplyr)

# Chemin vers votre fichier VCF
vcf_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/04-variants_calling/filtered_variants_subsamples.vcf"

# Lire le fichier VCF
vcf <- readVcf(vcf_path, "hg19")

# Afficher un résumé
summary(vcf)
# Afficher les premières lignes des variants
head(rowRanges(vcf))



#####Computing number of variant per super-scaffolds(all variants)#####
#######################################################################

###Infos on Quality###
#Variant quality values, represented by the QUAL field in a VCF (Variant Call Format) file, 
#provide a quantitative measure of confidence in the existence of the variant called at a specific position in the genome.  
#Variant quality is a measure of confidence that the variant (SNP, insertion, deletion, etc.) actually exists 
#in the sample analyzed, compared with a reference sequence. It is often expressed in phred score, which is a logarithmic scale.

# Extraire les séquences des noms directement depuis l'objet VCF
scaffolds <- seqnames(rowRanges(vcf))

# Compter le nombre de variants par super-scaffold
scaffold_counts <- table(scaffolds)

#creating data frame
df <- as.data.frame(scaffold_counts)
names(df) <- c("SuperScaffold", "VariantCount")

ggplot(df, aes(x = reorder(SuperScaffold, -VariantCount), y = VariantCount)) +
  geom_bar(stat = "identity", fill = "#AC6A9F") +
  theme_minimal() +
  labs(title = "Nombre de variants par super-scaffold", x = "Super-Scaffold", y = "Nombre de Variants") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))




#####Computing number of variant per super-scaffolds(over 50 quality)#####
##########################################################################


# Filtrer les variants avec une qualité inférieure à 50
filt_vcf <- subset(vcf, QUAL>50)

# Extraire les séquences des noms directement depuis l'objet VCF filtré
scaffolds_filt <- seqnames(rowRanges(filt_vcf))

# Compter le nombre de variants par super-scaffold après filtrage
scaffold_counts_filt <- table(scaffolds_filt)
head(scaffold_counts_filt)

# Créer un dataframe pour la visualisation des données filtrées
df_filt <- as.data.frame(scaffold_counts_filt)
names(df_filt) <- c("SuperScaffold", "VariantCount")

# Visualiser le nombre de variants par super-scaffold avec une qualité > 50
ggplot(df_filt, aes(x = reorder(SuperScaffold, -VariantCount), y = VariantCount)) +
  geom_bar(stat = "identity", fill = "#18A777") +
  theme_minimal() +
  labs(title = "Nombre de variants par super-scaffold (Qualité > 50)", 
       x = "Super-Scaffold", 
       y = "Nombre de Variants") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



#####Computing number of variant per super-scaffolds(over 100 quality)#####
##########################################################################


# Filtrer les variants avec une qualité inférieure à 50
filt_vcf <- subset(vcf, QUAL>100)

# Extraire les séquences des noms directement depuis l'objet VCF filtré
scaffolds_filt <- seqnames(rowRanges(filt_vcf))

# Compter le nombre de variants par super-scaffold après filtrage
scaffold_counts_filt <- table(scaffolds_filt)
head(scaffold_counts_filt)

# Créer un dataframe pour la visualisation des données filtrées
df_filt <- as.data.frame(scaffold_counts_filt)
names(df_filt) <- c("SuperScaffold", "VariantCount")

# Visualiser le nombre de variants par super-scaffold avec une qualité > 50
ggplot(df_filt, aes(x = reorder(SuperScaffold, -VariantCount), y = VariantCount)) +
  geom_bar(stat = "identity", fill = "#AE4C3C") +
  theme_minimal() +
  labs(title = "Nombre de variants par super-scaffold (Qualité > 100)", 
       x = "Super-Scaffold", 
       y = "Nombre de Variants") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))








