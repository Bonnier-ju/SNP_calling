
############################# Baypass results analysis ###############################
######################################################################################


################ Estimate the population covariance matrix ######################

library(ggplot2)
library(corrplot)
library(knitr)      
library(tidyverse)  
library(readxl)     
library(xtable)
library(reshape2)
library(kableExtra)
library(here)
library(parallel)
library(magrittr)
library(janitor)
library(readr)

# Upload the estimated Omega matrix
# The omega matrix quantifies how random genetic variations are correlated between different populations, 
# providing information on the genetic relationships and demographic structure of the populations studied.

omega_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.1-BayPass/Var_temp_pluvio/_mat_omega.out"
omega <- as.matrix(read.table(omega_file, header = F))
print(head(omega))

pop_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.1-BayPass/popfile.txt"
pop_names <- read_table(pop_file, col_names = FALSE)
pop_names <- as.character(pop_names$X1)

dimnames(omega) <- list(pop_names, pop_names)

# Visualization of the matrix
# as a correlation plot
cor_mat=cov2cor(omega)
corrplot(cor_mat,method="color",mar=c(2,1,2,2)+0.1, main=expression("Correlation map based on"~hat(Omega)))

# as a heatmap and hierarchical clustering tree (using the average agglomeration method)
hclust_ave <- function(x) hclust(x, method="average")
heatmap(1-cor_mat,hclustfun = hclust_ave,
        main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))


#################### Analysis of covariate model results of Baypass #######################

library(dplyr)


# extract SNP names
file.path.SNP_info <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/07-outlier_analysis/07.1-Bayescan/SNP_info.txt"
SNP_info=read.table(file.path.SNP_info)

# read result file
file.path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.1-BayPass/Run_2/SNP_run2_summary_betai_reg.out"
summary_data <- read.table(file.path, header = TRUE, sep = "")
head(summary_data)

colnames(summary_data)[1] <- "Env_var"

#split result file by covariables
data_covariable <- dplyr::filter(summary_data, Env_var == 14)

data_covariable$Env_var <- "BIO2 = Mean Diurnal Range"

data_covariable <- cbind(SNP_info, data_covariable)

names(data_covariable)[1:3] <- c("Scaffold", "Position", "SNP_name")
data_covariable[[5]] <- NULL
head(data_covariable)

# Extraire les lignes où BF.db. est supérieur à 20

filtered_data <- subset(data_covariable, BF.dB. > 20)

#filtered_data_BIO2 <- subset(data_covariable, BF.dB. > 20)

# Afficher les premières lignes pour vérifier
head(filtered_data)

# Sauvegarder le dataframe au format texte avec des tabulations comme séparateur
save_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.1-BayPass/Var_temp_pluvio/results_covariables/result_BIO12_Annual_Precipitation.txt"
write.table(filtered_data, file = save_path, sep = "\t", row.names = FALSE, quote = FALSE)


