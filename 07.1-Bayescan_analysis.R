
############################## Bayescan analysis ####################################
#####################################################################################

#r script created from : https://github.com/laurabenestan/Bayescan

library(ggplot2)

#Open the bayescan output file with the "_fst.txt" extension.

file.path.bayscan <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/07-outlier_analysis/07.1-Bayescan/bayescan_full_output_fst.txt"
file.path.SNP_info <- "C:/Users/bonni/OneDrive/Université/Thèse/Dicorynia/Article - Population Genomics/Analysis/07-outlier_analysis/07.1-Bayescan/SNP_info.txt"

bayescan=read.table(file.path.bayscan) 
SNP_info=read.table(file.path.SNP_info)

#The column of the bayescan object are prob, log10(P0), and qval are related to the test of local adaptation 
#considering the logarithm of the posterior odds - log10(PO) - and the q-value for the model with selection. 
#The fifth column gives the size of the locus-specific effect (alpha parameter). 
#The last one provides the locus-specific FST averaged over all populations.

#Combine SNP infos and bayescan results
bayescan = cbind(SNP_info, bayescan)

colnames(bayescan)[1:3] = c("Scaffold", "Position", "Name")

#Add a column for the type of selection grouping based on a Q-VALUE < 0.05 (you can also choose a Q-VALUE < 0.01).

bayescan$Selection <- ifelse(bayescan$alpha>=0&bayescan$qval<=0.01,"diversifying",ifelse(bayescan$alpha>=0&bayescan$qval>0.01,"neutral","balancing")) 
bayescan$Selection<- factor(bayescan$Selection)
levels(bayescan$Selection) 

#An ALPHA greater than zero suggests that allele frequency is increasing as a result of selection, 
#while an ALPHA less than zero indicates that selection favors the maintenance of several different alleles 
#at the same locus (balancing selection).
#A low Q_VALUE (e.g. < 0.05 or < 0.01) indicates that the results are statistically significant and therefore less 
#likely to be due to chance. Using thresholds such as 0.05 or 0.01 for Q_VALUE helps to control the rate of false 
#positives in the results, enabling more reliable interpretation of the data.

#Save the results of the SNPs potentially under positive (divergent) and balancing selection (qvalue < 0.05).
positive <- bayescan[bayescan$Selection=="diversifying",]
neutral <- bayescan[bayescan$Selection=="neutral",] 
balancing <- bayescan[bayescan$Selection=="balancing",]  

#Check the number of SNPs belonging to each category.
xtabs(data=bayescan, ~Selection)

#Write the results of the SNPs potentially under selection (qvalue < 0.05).
write.table(neutral, "neutral.txt", row.names=F, quote=F)  
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F)
