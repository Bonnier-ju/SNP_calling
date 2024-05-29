
########################## Redundancy analysis ################################
###############################################################################

#Loading package 
library(vcfR)
library(vegan)
library(dplyr)
library(RColorBrewer)
library(qvalue)
library(robust)
library(corrplot)
library(ggplot2)
library(car)

############################## Importing files and pre-processing per sites #######################################
###################################################################################################################

##### Genotype files #####

# Importing files
vcf_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/VCF_files/10per_04LD_pruned.vcf"
pop_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.4- RDA/Pop_id.csv"

vcf <- read.vcfR(vcf_file)
pop_data <- read.csv(pop_file)

# Extract genetic information
genotype_matrix <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

# Associate individuals with their populations
indiv_names <- colnames(genotype_matrix)
populations <- pop_data[match(indiv_names, pop_data$Sample_ID), "Population"]

# Prepare the data for aggregate function
genotype_data <- as.data.frame(t(genotype_matrix))
genotype_data$Population <- populations

# Use aggregate to calculate allele frequencies
AllFreq <- aggregate(genotype_data[, -ncol(genotype_data)], by = list(genotype_data$Population), FUN = function(x) mean(x, na.rm = TRUE) / 2)

######### The final genetic matrix included allele frequencies for 11 populations (rows) and loci (columns) ##########

##### Environmental file #####

# Loading environmental variables
env_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.4- RDA/Env_var_pop.csv"
env_data <- read.csv(env_file, row.names = 1)

# we standardized the climatic variables to avoid discrepancy in mean and standard deviation among variables and ensure that 
# the variable units were comparable (e.g., precipitation in mm and temperature in °C). 
# Note that RDA analyses can include qualitative predictors such as soil type using binary ("dummy variable") coding. 

## Standardization of the variables
Env_data <- scale(env_data, center=TRUE, scale=TRUE)

## Recovering scaling coefficients
scale_env <- attr(Env_data, 'scaled:scale')
center_env <- attr(Env_data, 'scaled:center')

## Climatic table
Env <- as.data.frame(Env_data)

# Loading geographic variables
geo_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.4- RDA/Geo_pop.csv"
geo_data <- read.csv(geo_file, row.names = 1)

# Standardization of the geographic variables
Geo_data <- scale(geo_data, center=TRUE, scale=TRUE)
scale_geo <- attr(Geo_data, 'scaled:scale')
center_geo <- attr(Geo_data, 'scaled:center')

# Convert to data frame and combine all variables
Geo <- as.data.frame(Geo_data)

Predictors <- cbind(Env, Geo)




############################## Importing files and pre-processing per individuals #################################
###################################################################################################################

# Importing files
vcf_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/VCF_files/10per_04LD_pruned.vcf"
pop_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.4- RDA/Pop_id.csv"

vcf <- read.vcfR(vcf_file)
pop_data <- read.csv(pop_file)

# Extract genetic information
genotype_matrix <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

# Association des individus avec leurs populations
indiv_names <- colnames(genotype_matrix)
populations <- pop_data[match(indiv_names, pop_data$Sample_ID), "Population"]

# Préparation des données pour calculer les fréquences alléliques par individu
genotype_data <- as.data.frame(t(genotype_matrix))
genotype_data$Individual <- indiv_names
genotype_data$Population <- populations

# Calcul des fréquences alléliques par individu
# Supposons que les valeurs des génotypes sont codées comme 0, 1, ou 2 (pour les génotypes homozygotes et hétérozygotes)
# La fréquence allélique pour un individu est la moyenne des valeurs divisée par 2 (pour normaliser entre 0 et 1)
AllFreqInd <- as.data.frame(lapply(genotype_data[, -c(ncol(genotype_data), ncol(genotype_data) - 1)], function(x) x / 2))

# Ajout des informations sur les individus et les populations
AllFreqInd$Individual <- genotype_data$Individual
AllFreqInd$Population <- genotype_data$Population


# Loading all variables per individuals
var_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.4- RDA/ALL_env_var_ids.csv"
Variables <- read.csv(var_file, row.names = 1)
Variables <- Variables[, -c(1)]

## Standardization of the variables
All_var <- scale(Variables, center=TRUE, scale=TRUE)

## Recovering scaling coefficients
scale_env <- attr(All_var, 'scaled:scale')
center_env <- attr(All_var, 'scaled:center')

## Climatic table
All_env_var <- as.data.frame(All_var)

##### Goegraphic variables #####

# Loading geographic variables
geo_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.4- RDA/Geo_ids.csv"
geo_data <- read.csv(geo_file, row.names = 1)

# Standardization of the geographic variables
Geo_data <- scale(geo_data, center=TRUE, scale=TRUE)
scale_geo <- attr(Geo_data, 'scaled:scale')
center_geo <- attr(Geo_data, 'scaled:center')

# Convert to data frame and combine all variables
Geo <- as.data.frame(Geo_data)

Predictors <- cbind(All_env_var, Geo)


############################## Environmental variables processing ###############################################
#################################################################################################################

############################## Variable selection: forward model building procedure #############################


# Forward selection starts from a "null" model where the response is explained only by an intercept. 
# Variables are then added to the model one by one to try to reach the amount of variance explained 
# by a "full" model (i.e., model including all the explanatory variables), while limiting the amount 
# of redundancy among included variables.    

## Null model
RDA0 <- rda(AllFreqInd[, -1] ~ 1, All_env_var)

## Full model
RDAfull <- rda(AllFreqInd[, -1] ~ 
                 Pedologie + 
                 Mean.annual.temperature + 
                 Mean.Diurnal.Range + 
                 Isothermality + 
                 Temperature.Seasonality + 
                 Max.Temperature.of.Warmest.Month + 
                 Min.Temperature.of.Coldest.Month + 
                 Temperature.Annual.Range + 
                 Mean.Temperature.of.Wettest.Quarter + 
                 Mean.Temperature.of.Driest.Quarter + 
                 Mean.Temperature.of.Warmest.Quarter + 
                 Mean.Temperature.of.Coldest.Quarter + 
                 Annual.Precipitation + 
                 Precipitation.of.Wettest.Month + 
                 Precipitation.of.Driest.Month + 
                 Precipitation.Seasonality..Coefficient.of.Variation. + 
                 annualPET + 
                 aridityIndex + 
                 climaticMoisture + 
                 continentality + 
                 maxTempColdestMonth + 
                 minTempWarmestMonth + 
                 PETColdestQuarter + 
                 PETDriestQuarter + 
                 PETseasonality + 
                 PETWarmestQuarter + 
                 PETWettestQuarter + 
                 thermInd + 
                 tri + 
                 topoWet, 
               All_env_var)

# To conduct the selection procedure we used the ordiR2step function of the package vegan 
# and the following stopping criteria: variable significance of p < 0.01 using 1000 permutations, 
# and the adjusted R2 of the global model.  


## Stepwise procedure with ordiR2step function
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)

mod$anova



########################### Processing PCA on genotype #################################
########################################################################################


# Perform PCA on the aggregated allele frequencies (neutral genetic markers)
pca_result <- rda(AllFreq[, -1], scale = TRUE) # PCA in vegan uses the rda() call without any predictors

# Screeplot of the PCA eigenvalues
screeplot(pca_result, type = "barplot", npcs = 10, main = "PCA Eigenvalues")

# Extract the first three principal components
PCs <- scores(pca_result, choices = c(1:3), display = "sites", scaling = 0)
PopStruct <- data.frame(Population = AllFreq[, 1], PCs)
colnames(PopStruct) <- c("Population", "PC1", "PC2", "PC3")

ggplot(PopStruct, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Population Structure PCA", x = "PC1", y = "PC2")

# Combine environmental, geographic, and PCA axes data
# Ensure Population column is removed before combining
Predictors <- cbind(Env, Geo, PopStruct[, -1])



############################### Performing partial RDA analysis ##########################
##########################################################################################

selected_predictors <- Predictors[, c("PC1", "PC2", "PC3", "long", "lat", "MoistureIndex", "annualPET", "PETWettestQuarter", "Continentality")]

# Full model
pRDAfull <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + long + lat + MoistureIndex + annualPET + PETWettestQuarter + Continentality, data = selected_predictors)
full_adj_r2 <- RsquareAdj(pRDAfull)
full_anova <- anova(pRDAfull)

vif(pRDAfull)

# Pure env model
pRDAclim <- rda(AllFreq[, -1] ~ MoistureIndex + annualPET + PETWettestQuarter + Continentality + Condition(long + lat + PC1 + PC2 + PC3), data = selected_predictors)
clim_adj_r2 <- RsquareAdj(pRDAclim)
clim_anova <- anova(pRDAclim)

# Pure neutral population structure model
pRDAstruct <- rda(AllFreq[, -1] ~ PC1 + PC2 + PC3 + Condition(long + lat + MoistureIndex + annualPET + PETWettestQuarter + Continentality), data = selected_predictors)
struct_adj_r2 <- RsquareAdj(pRDAstruct)
struct_anova <- anova(pRDAstruct)

# Pure geography model
pRDAgeog <- rda(AllFreq[, -1] ~ long + lat + Condition(MoistureIndex + annualPET + PETWettestQuarter + Continentality + PC1 + PC2 + PC3), data = selected_predictors)
geog_adj_r2 <- RsquareAdj(pRDAgeog)
geog_anova <- anova(pRDAgeog)

# Summary of results 
list(
  Full_Model = list(Adj_R2 = full_adj_r2, ANOVA = full_anova),
  Climate_Model = list(Adj_R2 = clim_adj_r2, ANOVA = clim_anova),
  Structure_Model = list(Adj_R2 = struct_adj_r2, ANOVA = struct_anova),
  Geography_Model = list(Adj_R2 = geog_adj_r2, ANOVA = geog_anova)
)

# View correlations among selected variables
corrplot(cor(selected_predictors), type = "upper")

#**Notes on interpretation and best practices:** 
# In this case, the largest proportion of genetic variance could not be uniquely attributed to any of the three sets 
# of predictors, a common occurrence given the ubiquitous nature of spatial autocorrelation in environmental and genetic data 
# sets. This confounded effect reflects a high degree of collinearity among explanatory variables. 
# This is critical information given that most landscape genomic studies look for correlation between climatic and genetic variation (i.e., GEA) 
# and either assume no collinearity or, on the contrary, totally remove this commonly explained variation. 
# In the first case, GEA detections could potentially be subject to high false positive rates, while in the latter case detections 
# might show high false negative rates. Selecting an appropriate approach to account for demographic history and geographic distance 
# is of major importance when searching for selection in the genome. Variance partitioning can be a useful step to explore 
# the (statistical) association among available descriptors, to better understand the covariation of environmental and genetic gradients, 
# and to determine how much overall genetic variation is shaped by environmental, geographic, and demographic factors before conducting further 
# landscape genomics study.




############################ Genotype-Environment Associations: identifying loci under selection ##############################
###############################################################################################################################

#### Function to conduct a RDA based genome scan
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}


# The first step was to run an RDA model on the allele frequency matrix using the retained climatic factors 
# as explanatory variables and the first three PCs as conditioning variables to account for neutral population structure.  
RDA_env <- rda(AllFreq[, -1] ~ MoistureIndex + annualPET + PETWettestQuarter + Continentality + Pedologie + BIO15...Precipitation.Seasonality..Coefficient.of.Variation. + Condition(PC1 + PC2 + PC3), data = Predictors)

# We then had to choose a number of RDA axes to include when conducting the genome scan.  
screeplot(RDA_env, main="Eigenvalues of constrained axes")

## Running the function with K = 3
rdadapt_env<-rdadapt(RDA_env, 3)

plot(-log10(rdadapt_env$p.values), main = "Manhattan Plot des p-values", xlab = "Position SNP", ylab = "-log10(p-value)")
abline(h = -log10(thres_env), col = "red")

# One critical step when conducting a genome scan is to set a pertinent p-value threshold to identify the outlier loci. 
# Here, we used a Bonferroni correction to account for multiple testing.  

# the rdadapt function returns both p-values and q-values, which means it is possible to use a FDR (False Discovery Rate) 
# approach instead of a p-value threshold to identify outliers

## P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values)

## Identifying the loci that are below the p-value threshold
outliers_pval <- data.frame(
  Loci = colnames(AllFreq)[which(rdadapt_env$p.values < thres_env)], 
  p.value = rdadapt_env$p.values[which(rdadapt_env$p.values < thres_env)]
)

fdr_threshold <- 0.05/length(rdadapt_env$q.values)

## Identifying the loci that are below the FDR threshold
outliers_fdr <- data.frame(
  Loci = colnames(AllFreq)[which(rdadapt_env$q.values < fdr_threshold)], 
  q.value = rdadapt_env$q.values[which(rdadapt_env$q.values < fdr_threshold)]
)


# Once the outliers have been identified, it can be useful to visualize their distribution in comparison with 
# neutral loci using either an RDA biplot or a Manhattan plot.

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
#TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))



