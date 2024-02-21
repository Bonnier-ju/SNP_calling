# Load necessary libraries
library(tidyverse)
library(cowplot)

# Define the base path for the input files
base_path <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Bio-informatique analysis/05-filtering_variants/05.2-missing_SNP"

# Load population information from a CSV file
pop_info <- read_csv("Population_info.csv")

# Analyzing missing data by individual
# Load the missing data statistics for individuals
missing_indv_path <- file.path(base_path, "all.biallelic.snp.filtered.imiss")
missing_indv <- read_delim(missing_indv_path, delim = " ") %>%
  rename_all(~gsub(" ", "", .)) %>%
  mutate_at(vars(N_MISS, N_GENO, F_MISS), as.numeric) %>%
  mutate(Ind = gsub(".g.vcf", "", IID)) %>%
  mutate(Ind = gsub(" ", "", Ind)) %>%
  left_join(pop_info, by = c("Ind" = "SampleID"))

# Visualize the proportion of missing SNPs per individual with a histogram
g.het_imiss <- ggplot(missing_indv, aes(F_MISS, fill = Pop)) +
  geom_histogram(position = "dodge") +
  scale_x_log10() +
  xlab("Proportion of missing SNPs per individual") +
  ylab("Number of individuals")

# Analyzing PCA based on the filtered biallelic SNP data
pca_path <- file.path(base_path, "symcapture.all.biallelic.snp.filtered.eigenvec")
g.pca <- read_delim(pca_path, delim = " ", col_names = FALSE) %>%
  rename(Sample = V2, PCA1 = V3, PCA2 = V4) %>%
  mutate(Ind = gsub(".g.vcf", "", Sample)) %>%
  left_join(pop_info, by = c("Ind" = "SampleID")) %>%
  ggplot(aes(x = PCA1, y = PCA2, color = Pop)) + 
  geom_point(size=2) +
  stat_ellipse(level = 0.95, size = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

# Plot both the missing data histogram and PCA plot together
cowplot::plot_grid(g.het_imiss, g.pca, nrow = 2, rel_heights = c(1,2))

# Analyzing missing data by SNP
missing_snp_path <- file.path(base_path, "symcapture.all.biallelic.snp.filtered.lmiss")
missing_snp <- read_delim(missing_snp_path, delim = " ") %>%
  rename_all(~gsub(" ", "", .)) %>%
  mutate(F_MISS = as.numeric(F_MISS)) %>%
  ggplot(aes(F_MISS)) +
  geom_histogram(position = "dodge") +
  xlab("Proportion of missing data per SNP") +
  ylab("SNP count") +
  scale_y_sqrt() +
  scale_x_sqrt()

print(missing_snp)

