
################### Basics statistics on population ############################
################################################################################

library(vcfR)
library(hierfstat)
library(adegenet)


# Files path 
vcf_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/VCF_files/sub_200k_SNP_stgeorge.vcf"
#pop_file_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/VCF_files/Pop_file_East_Inland_West_StG.csv"

# Read data files 
vcf_data <- read.vcfR(vcf_path)
#pop_file <- read.csv(pop_file_path, header = T)

# Converting VCF to genind 
genind_obj <- vcfR2genind(vcf_data)
genind_obj

# Creating vecteur for pop name
groups <- c("Inland", "Inland", "East", "East", "Inland", "East", "East", "East", "Inland", "East", 
                 "East", "St_georges", "East", "East", "East", "East", "East", "Inland", "St_georges", "Inland", 
                 "East", "East", "West", "East", "East", "West", "West", "Inland", "East", "West", "St_georges", 
                 "East", "St_georges", "Inland", "Inland", "West", "Inland", "Inland", "East", "Inland", "St_georges", 
                 "Inland", "Inland", "East", "East", "Inland", "Inland", "East", "Inland", "West", "West", "Inland", 
                 "Inland", "East", "Inland", "Inland", "Inland", "West", "Inland", "Inland", "Inland", "Inland", "East", 
                 "West", "Inland", "West", "Inland", "West", "West", "West", "Inland", "Inland", "West", "East", "East", 
                 "West", "St_georges", "St_georges", "Inland", "West", "Inland", "Inland", "East", "Inland", "Inland", 
                 "Inland", "St_georges")

groups <- c("St_georges", "St_georges", "St_georges", "St_georges", "St_georges", "St_georges", "St_georges", "St_georges")

# Check if same number of individuals and pop name
if (length(groups) != nInd(genind_obj)) {
  stop("Le nombre de populations ne correspond pas au nombre d'individus dans l'objet genind.")
}

# Add pop to genind
pop(genind_obj) <- factor(groups)

table(pop(genind_obj))

#### Basic stats by hierfstats ####
#Estimates individual counts, allelic frequencies, observed heterozygosities and genetic diversities per locus and population. 
#Also Estimates mean observed heterozygosities, mean gene diversities within population Hs, Gene diversities overall Ht and corrected Htp, and Dst, Dstp. 
#Finally, estimates Fst and Fstp as well as Fis following Nei (1987) per locus and overall loci

stats <- basic.stats(genind_obj)

# Saving results in differents data-frame
stats_df <- do.call(rbind, stats)
str(stats)

# Convert the numeric matrices to data frames
n_ind_samp_df <- as.data.frame(stats$n.ind.samp)
Ho_df <- as.data.frame(stats$Ho)
Hs_df <- as.data.frame(stats$Hs)
Fis_df <- as.data.frame(stats$Fis)

perloc_df <- stats$perloc

# Convert the 'overall' named vector to a data frame
overall_df <- as.data.frame(t(stats$overall))
names(overall_df) <- names(stats$overall)

# Processing pop.freq tables into a single data frame
pop_freq_df <- do.call(rbind, lapply(names(stats$pop.freq), function(marker) {
  freq_table <- as.data.frame.matrix(stats$pop.freq[[marker]])
  freq_table$Marker = marker  # Add marker identifier
  return(freq_table)
}))

# Write data frames to CSV files
write.csv(n_ind_samp_df, "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.3-basic_stats/n_ind_samp.csv", row.names = FALSE)
write.csv(Ho_df, "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.3-basic_stats/Ho.csv", row.names = FALSE)
write.csv(Hs_df, "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.3-basic_stats/Hs.csv", row.names = FALSE)
write.csv(Fis_df, "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.3-basic_stats/Fis.csv", row.names = FALSE)
write.csv(perloc_df, "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.3-basic_stats/perloc.csv", row.names = FALSE)
write.csv(overall_df, "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.3-basic_stats/overall.csv", row.names = FALSE)
write.csv(pop_freq_df, "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/06-population_structure/06.3-basic_stats/pop_freq.csv", row.names = FALSE)

#Summary per groups
summary(Ho_df)
summary(Hs_df)
summary(Fis_df)
summary(overall_df)
summary(pop_freq_df)


###################### Test of Hardy--Weinberg Equilibrium #####################
################################################################################
#This function tests, for a series of loci, the hypothesis that genotype frequencies follow the Hardy--Weinberg equilibrium.

hw.test(genind_obj, B = 1000)




##################### Allelic richness by hierfstats ###########################
################################################################################
#Estimates allelic richness, the rarefied allelic counts, per locus and population

A_richness <- allelic.richness(genind_obj)

str(A_richness)

Ar_df <- as.data.frame(A_richness$Ar)

summary(Ar_df)





############################# DAPC ############################################
###############################################################################

library(adegenet)
library(ade4)

#color package 
library(RColorBrewer)
library(viridis)


vcf_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/VCF_files/sub_200k_SNP.vcf"

vcf_data <- read.vcfR(vcf_path)

genind_obj <- vcfR2genind(vcf_data)

# Creating a vector for pop name
groups <- c("Cacao","Cacao","MC_87","MC_88","Cacao","MC_88","MC_88","MC_87","Regina","MC_88","Foret_Regina_St_Georges","St_georges",
            "Saut_Lavilette","Saut_Lavilette","Saut_Lavilette","MC_87","Saut_Lavilette","Nouragues_Inselberg","St_georges","Regina",
            "Foret_Regina_St_Georges","Foret_Regina_St_Georges","Apatou","Foret_Regina_St_Georges","Saut_Lavilette","Apatou","Apatou",
            "Regina","MC_88","Acarouany","St_georges","MC_88","St_georges","Nouragues_Inselberg","Piste_St_Elie","Acarouany","Regina",
            "Nouragues_Inselberg","MC_88","Piste_St_Elie","St_georges","Regina","Nouragues_Inselberg","Foret_Regina_St_Georges","Saut_Lavilette",
            "Nouragues_Inselberg","Cacao","MC_88","Piste_St_Elie","Acarouany","Apatou","Nouragues_Inselberg","Nouragues_Inselberg",
            "Foret_Regina_St_Georges","Piste_St_Elie","Piste_St_Elie","Nouragues_Inselberg","Apatou","Cacao","Cacao","Nouragues_Inselberg",
            "Nouragues_Inselberg","MC_87","Apatou","Piste_St_Elie","Acarouany","Nouragues_Inselberg","Acarouany","Acarouany","Acarouany",
            "Nouragues_Inselberg","Nouragues_Inselberg","Acarouany","Saut_Lavilette","Saut_Lavilette","Apatou","St_georges","St_georges",
            "Regina","Apatou","Nouragues_Inselberg","Regina","Foret_Regina_St_Georges","Piste_St_Elie","Nouragues_Inselberg","Regina","St_georges")

# Check that the length of the population vector corresponds to the number of individuals in genind_obj
if (length(groups) != nInd(genind_obj)) {
  stop("Le nombre de populations ne correspond pas au nombre d'individus dans l'objet genind.")
}

# Assign pop to genind
pop(genind_obj) <- factor(groups)

table(pop(genind_obj))

# Find the optimal number of clusters without population a priori
# Find the break int he curve for number of retained PCA 
# Use the lowest value of number of cluster for BIC values 
grp <- find.clusters(genind_obj, max.n.clust = 30)

# Perform Principal Component Analysis (PCA) to reduce data dimensionality
# and avoid statistical overdetermination
dapc_obj <- dapc(genind_obj, var.contrib = TRUE, n.pca = 50, n.da = length(unique(groups)) - 1)

# Link inds with pop
dapc_obj$grp <- factor(groups)

# Create a color vector where names correspond to groups
specific_colors <- c("Acarouany" = "#D02090", "Apatou" = "red", "Cacao" = "green", 
                     "Foret_Regina_St_Georges" = "#00FFFF", "MC_87" = "chartreuse", "MC_88" = "#008B8B",
                     "Nouragues_Inselberg" = "deepskyblue", "Piste_St_Elie" = "chartreuse4", 
                     "Regina" = "dodgerblue2", "Saut_Lavilette" = "royalblue3", "St_georges" = "#5CACEE")

myCol <- c("darkblue","purple","green","orange","red","blue")

# DAPC
scatter(dapc_obj, bg="white", pch=20, cell=4, cstar=2, col = specific_colors[dapc_obj$grp], solid=1,cex=1,clab=0.5, leg=TRUE, scree.pca=TRUE, scree.da=TRUE,
        posi.pca="bottomleft", posi.da="bottomright", ratio.pca=0.16, ratio.da=0.16, cleg = .6)

# Calculation of the percentage of variance explained by the main axes
eigenvalues <- dapc_obj$eig
percent_variance <- eigenvalues / sum(eigenvalues) * 100




######################################### FST computing #########################################################
#################################################################################################################
# Install the necessary Bioconductor packages
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("SNPRelate")

install.packages("dartRverse")
install.packages("dartR")

library(dartRverse)
library(dartR)
library(vcfR)
library(adegenet)


vcf_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/VCF_files/ld_04_pruned.vcf"
csv_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/Pop_files/Pop_file_11sites.csv"


data <- read.vcfR(vcf_path)

data <- vcfR2genind(data)

pop_info <- read.csv(csv_file, header = T)

populations <- pop_info[, "population"]

genind_data@pop <- factor(populations)

print(genind_data)

data <- gi2gl(data, parallel = FALSE, verbose = NULL)
