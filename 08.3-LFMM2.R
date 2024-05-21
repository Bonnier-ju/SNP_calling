
################################# Analysis with LFMM2 ######################################
# The latent factor mixed model (LFMM) is a multivariate mixed regression model that estimates simultaneously the effects of environmental 
# variables (fixed effects) and unobserved confounders called latent factors. The latent factors are computed both from the genomes 
# and from their environment. They are not representing neutral population structure (i.e. they have less direct interpretations than in 
# ancestry estimation methods). Instead, they can be interpreted as the best estimates of the confounding effects of neutral population 
#structure, leading to environmental effect size estimates with minimal bias.



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

if(!requireNamespace("qvalue", quietly = TRUE)) {  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("qvalue")
}

library(LEA)
library(vcfR)
library(qvalue)

###################### Pre-processsing of SNP data #############################

vcf_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/VCF_files/10per_04LD_pruned.vcf"
lfmm_output_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/VCF_files/10per_04LD_pruned.lfmm"

vcf2lfmm(input.file = vcf_file, output.file = lfmm_output_file, force = TRUE)
genotype_matrix <- read.lfmm(lfmm_output_file)


###################### Pre-processsing of env data #############################
# We attribute climatic values to each genotype, i.e. genotypes from the same populations will have the same climatic values.

env_csv_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.3-LFMM/Var_env_by_inds.csv"
env_file <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/08-EAA/08.3-LFMM/env_file.env"

env_data <- read.csv(env_csv_file, header = TRUE)[, -c(1, 2)]
head(env_data)
write.table(env_data, env_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
env_data <- read.table(env_file, sep = "\t", header = TRUE)


nrow(genotype_matrix)
nrow(env_data)

############################ Running LFMM2 ########################################

# LFMM requires an estimate of the number of populations in the data (K)
# Following our previous analysis on population structure best K is 3 
K <- 3

# Run LFMM2
# The lfmm2 function estimates latent factors based on an exact least-squares approach.
# The resulting object can be used by the function lfmm2.test to identify genetic polymorphisms 
# exhibiting association with ecological gradients or phenotypes, while correcting for unobserved confounders. 
mod_lfmm2 <- lfmm2(input = genotype_matrix, env = env_data, K = K)

# uncalibrated p-values
# With the lfmm2.test function, we can obtain a vector of p-values for associations between loci and 
# environmental variables adjusted for latent factors calculated by lfmm2.
lfmm2.test(object = mod_lfmm2, input = genotype_matrix, env = env_data, full = FALSE, genomic.control = FALSE)$pvalues %>%
  hist(col = "orange", 
       main="Histogram of non-calibrated p-values",
       xlab="p-values")

# p-values calibrated
lfmm2.test(object = mod_lfmm2, input = genotype_matrix, env = env_data, full = FALSE, genomic.control = TRUE)$pvalues %>%
  hist(col = "orange", 
       main="Histogram of calibrated p-values",
       xlab="p-values")

# Looking for the ideal p-value distribution by using calibrated p-value => flat with a peak at 0

test_lfmm2 <- lfmm2.test(object = mod_lfmm2, input = genotype_matrix, env = env_data, full = FALSE, genomic.control = TRUE)
pv_lfmm2 <- test_lfmm2$pvalues

# Renamin col of p-value file with true SNP names
vcf <- read.vcfR(vcf_file)
snp_names <- vcf@fix[, "ID"]  
head(snp_names)
colnames(pv_lfmm2) <- snp_names

# We convert the adjusted p-values to q-values. q-values provide a measure of each SNP’s significance, 
# automatically taking into account the fact that thousands are simultaneously being tested. 
# We can then use an FDR (False Discovery Rate) threshold to control the number of false positive detections 
# (given that our p-value distribution is “well-behaved”).

fdr_level <- 0.01

# Extract p-value for one variable
annual_pet_pvalues <- pv_lfmm2["annualPET", ]

# Apply FDR
qv_lfmm2 <- qvalue::qvalue(annual_pet_pvalues, fdr.level = fdr_level)

# Identify candidat SNP
significant_snps <- which(qv_lfmm2$qvalues < fdr_level)
# or 
# candidates <- which(qv_lfmm2$significant)

# Create à data frame
significant_snps_df <- data.frame(SNP = significant_snps, 
                                  p_value = valid_pvalues[significant_snps], 
                                  q_value = qv_lfmm2$qvalues[significant_snps])

head(significant_snps_df)


