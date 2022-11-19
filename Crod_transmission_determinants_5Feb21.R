library(ggplot2)
library(ggforce)
library(plyr)
library("dplyr")
library(tidyr)
library("beeswarm")
library("viridis")
library("fitdistrplus")
library("hrbrthemes")
library("ggtree")
library(cowplot)
library(tidyverse)
library("ggpubr")
library("ActuDistns")
library("vcd")
library("reshape2")
library(grid)
library(gridExtra)
library("svglite")

setwd("~/Projects/C._rodentium/gvcf_reanalysis/Final_analysis/")

dists <- read.csv("Citrobacter_distances_table.csv")
alleles <- read.csv("All_snps_summary_long_clean.csv")

#Assign query sites: sites with Zsnps or at least 6 hSNPs
#Remove the noisy sites
Good_sites <- c(2438536,3996563,1285412,2835878,793586,1039375,1229417,17805,1749372,1896011,4987380,492511,2372393,175076,2194913,4623423,645840,3040456,2584041,1417667,123360,890802,937620,2172181,408983,3638848,1531388,2169284,5081862,1086564,4298690,1415559,4854857,1032727,1265896,3928225)
#Sites deemed to be noise
Noisy <- c(123360,175076,408983,645840,890802,1086564,1417667,2169284,2172181,2194913,2584041,3040456,4623423,5081862)
#Sites with fixed SNPs
Zsnps <- c(2372393 , 492511 , 1896011 , 4987380 , 17805 , 1749372 , 1229417 , 1039375 , 793586 , 1285412 , 2835878 , 3996563 , 2438536)

#Store variable H SNP sites to a vector
hSNPs <- Good_sites[!Good_sites %in% Zsnps]

#Store variable H SNP sites to a vector without noisy sites
cleanhSNPS <- hSNPs[!hSNPs %in% Noisy]

# Identify hSNPS that are in ICC180
ICC_hSNPs <- alleles$Pos[which(alleles$ID == "ICC180" & alleles$Allele_Ratio > 0.025)]
 
##---------------------------------------------------------------------------
# Step 1 to characterize pairwise comparisons in the dists dataframe
##---------------------------------------------------------------------------
##---------------------------------------------------------------------------
# Add column for matched chains
for (i in 1:nrow(dists)){
  ch1 <- levels(dists$Ch1[i])[dists$Ch1[i]]
  ch2 <- levels(dists$Ch2[i])[dists$Ch2[i]]
  if (ch1 == ch2){
    dists$Same_chain[i] <- "Yes"
    
  }else {
    dists$Same_chain[i] <- "No"
  }
}

# Calculate the transmission time steps between each pair

for (i in 1:nrow(dists)){
  ch1 <- levels(dists$Ch1[i])[dists$Ch1[i]]
  ch2 <- levels(dists$Ch2[i])[dists$Ch2[i]]
  if (ch1 == ch2){
    dists$timestep[i] <- abs(dists$M1[i] - dists$M2[i])
    
  }else {
    dists$timestep[i] <- dists$M1[i] + dists$M2[i]
  }
}

##---------------------------------------------------------------------------
# Look at HSNP sites and ask the following question:
# 1) How many isolates have it as a hSNPs
# 2) How frequently is it shared between transmission pairs
# 3) how frequently is it shared between non transmission pairs
##---------------------------------------------------------------------------

hSites <- alleles %>% dplyr::select(Pos) %>% unique()


## Quantifies how many times a site is shared across transmission pairs, cluster pairs and chain pairs
for (i in 1:nrow(hSites)){
  hSites$frequency[i] <-  alleles %>% filter(Pos == hSites$Pos[i], Allele_Ratio > 0.025) %>% nrow()

  trans_hvmatches <- 0
  no_trans_hvmatches <- 0
  
  cluster_hvmatches <- 0
  no_cluster_hvmatches <- 0
  
  chain_hvmatches <- 0
  no_chain_hvmatches <- 0
  
  for(j in 1:nrow(dists)){
    ID1 <- levels(dists$ID1[j])[dists$ID1[j]]
    ID2 <- levels(dists$ID2[j])[dists$ID2[j]]
    index1 <- which(alleles$Pos == hSites$Pos[i] & alleles$ID == ID1 )
    index2 <- which(alleles$Pos == hSites$Pos[i] & alleles$ID == ID2 )
    
    if(dists$timestep[j] == 1){
      if(alleles$Allele_Ratio[index1] > 0.025 && alleles$Allele_Ratio[index2] > 0.025){
        trans_hvmatches <- trans_hvmatches + 1
      }
    } 
    if(dists$timestep[j] > 1){
      if(alleles$Allele_Ratio[index1] > 0.025 && alleles$Allele_Ratio[index2] > 0.025){
        no_trans_hvmatches <- no_trans_hvmatches + 1
      }
    }
    if(dists$timestep[j] < 6){
      if(alleles$Allele_Ratio[index1] > 0.025 && alleles$Allele_Ratio[index2] > 0.025){
        cluster_hvmatches <- cluster_hvmatches + 1
      }
    } 
    if(dists$timestep[j] > 5){
      if(alleles$Allele_Ratio[index1] > 0.025 && alleles$Allele_Ratio[index2] > 0.025){
        no_cluster_hvmatches <- no_cluster_hvmatches + 1
      }
    }
    if(dists$Same_chain[j] == "Yes"){
      if(alleles$Allele_Ratio[index1] > 0.025 && alleles$Allele_Ratio[index2] > 0.025){
        chain_hvmatches <- chain_hvmatches + 1
      }
    } 
    if(dists$Same_chain[j] == "No"){
      if(alleles$Allele_Ratio[index1] > 0.025 && alleles$Allele_Ratio[index2] > 0.025){
        no_chain_hvmatches <- no_chain_hvmatches + 1
      }
    }
  }
  
  hSites$trans_hvshared[i] <- trans_hvmatches
  hSites$Notrans_hvshared[i] <- no_trans_hvmatches
  
  hSites$cluster_hvshared[i] <- cluster_hvmatches
  hSites$Nocluster_hvshared[i] <- no_cluster_hvmatches
  
  hSites$chain_hvshared[i] <- chain_hvmatches
  hSites$Nochain_hvshared[i] <- no_chain_hvmatches
  
}


### Use the multiplicity rule to deter P(trans|shared) = P(trans & shared)/P(Shared)
#for(i in 1:nrow(hSites)){
#hSites$PTS[i] <- hSites$trans_shared[i]/(hSites$trans_shared[i] + hSites$Notrans_shared[i])
#hSites$PTCl[i] <- hSites$cluster_shared[i]/(hSites$cluster_shared[i] + hSites$Nocluster_shared[i])
#hSites$PTCh[i] <- hSites$chain_shared[i]/(hSites$chain_shared[i] + hSites$Nochain_shared[i])
#}

### Use Bayes Theorem to deter P(trans|shared) = P(shared|trans)*P(trans)/P(Shared)
ntrans <- nrow(dists %>% filter(timestep==1))
ncluster <- nrow(dists %>% filter(timestep < 6))
nchain <- nrow(dists %>% filter(Same_chain == "Yes"))

ptrans <- mean(dists$timestep == 1)
pcluster <- mean(dists$timestep < 6)
pchain <- mean(dists$Same_chain == "Yes")

for(i in 1:nrow(hSites)){
  hSites$transmission_predictiveness[i] <- (hSites$trans_hvshared[i]/ntrans)*ptrans/((hSites$trans_hvshared[i] + hSites$Notrans_hvshared[i])/nrow(dists))
  hSites$Cluster_predictiveness[i] <- (hSites$cluster_hvshared[i]/ncluster)*pcluster/((hSites$cluster_hvshared[i] + hSites$Nocluster_hvshared[i])/nrow(dists))
  hSites$Same_Chain_predictiveness[i] <- (hSites$chain_hvshared[i]/nchain)*pchain/((hSites$chain_hvshared[i] + hSites$Nochain_hvshared[i])/nrow(dists))
}

hSites$Designation <- "iSNV"
for(i in 1:nrow(hSites)){
  qs <- hSites$Pos[i]
  if(qs %in% Noisy){
    hSites$Designation[i] <- "Noisy"
  }
  if(qs %in% Zsnps){
    hSites$Designation[i] <- "SNV"
  }
}
##---------------------------------------------------------------------------
# Plot the predictiveness of sites
##---------------------------------------------------------------------------

PTS_hv_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = transmission_predictiveness,color = Designation)) + 
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  ylab("P(Transmission|Shared)") + xlab("") +
  geom_point() + theme(legend.position='none')

PTCl_hv_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = Cluster_predictiveness, color = Designation)) + 
  ylab("P(Cluster|Shared)") + xlab("Prevalence (%)") +
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  geom_point() + theme(legend.position='none')

PTCh_hv_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = Same_Chain_predictiveness, color = Designation, )) + 
  ylab("P(Same chain|Shared)") + xlab("") +
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  geom_point() + theme(legend.position="bottom") + theme(legend.title = element_blank())


x_blanker <- theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

legend <- as_ggplot(get_legend(PTCh_hv_plot))

PTCh_hv_plot <- PTCh_hv_plot + theme(legend.position='none')

predictive_sites_panel <- plot_grid(PTS_hv_plot,PTCl_hv_plot,PTCh_hv_plot, ncol =  3, labels = "AUTO")
predictive_sites_panel <- plot_grid(predictive_sites_panel,legend, nrow = 2, rel_heights = c(9,1))

Pred_tbl <- gather((hSites %>% filter(frequency > 2) %>% dplyr::select(Pos,frequency,transmission_predictiveness,Cluster_predictiveness, Same_Chain_predictiveness)), Metric, Predictiveness, transmission_predictiveness:Same_Chain_predictiveness, factor_key=TRUE)
Variant_sites <- hSites %>% filter(!Pos %in% Noisy) %>% filter(frequency != 0) %>% pull(Pos)
Variants_not_rare <- hSites %>% filter(!Pos %in% Noisy) %>% filter(frequency > 2) %>% pull(Pos)
write.csv(hSites, "Sites_predictiveness.csv", row.names = FALSE)

### Export image to pdf as Figure 1
save_plot("Sites_predictiveness.pdf",predictive_sites_panel , base_height = 4, base_width = 8)

##---------------------------------------------------------------------------
# For each comparison count how many times a new variant emerges
##---------------------------------------------------------------------------
trans_dists <- dists %>% filter(timestep == 1)

for (i in 1:nrow(trans_dists)) {
  ID1 <- levels(trans_dists$ID1[i])[trans_dists$ID1[i]]
  ID2 <- levels(trans_dists$ID2[i])[trans_dists$ID2[i]]
  new_sites <- 0
  
  for (qpos in Variant_sites ) {
    q1 <- which(alleles$Pos == qpos & alleles$ID == ID1 )
    q2 <- which(alleles$Pos == qpos & alleles$ID == ID2 )
    
    if(alleles$Allele_Ratio[q1] <= 0.025 & alleles$Allele_Ratio[q2] > 0.025 ){
      new_sites <- new_sites + 1
    }
  }
  trans_dists$New_variants[i] <- new_sites
}  

##---------------------------------------------------------------------------
# Compare the mutation rate and the substitution rate between N and W chains
##---------------------------------------------------------------------------

trans_dists$Treatment <- substr(trans_dists$Ch2,1,1)

SNPs_plot <- trans_dists  %>% ggplot(aes(x = Dist)) + 
  geom_histogram(bins = 3, aes(y=..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.5,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of SNPs per transmission") + ylab("Relative frequency") +
  theme(legend.position = "none")

SNPs_plot_Treatment <- trans_dists  %>% ggplot(aes(x = Dist)) + 
  geom_histogram(bins = 3, aes(y=2*..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.5,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of SNPs per transmission") + ylab("Relative frequency") +
  theme(legend.position = "none") + facet_wrap(~Treatment)

variants_plot <- trans_dists %>% ggplot(aes(x = New_variants)) + 
  geom_histogram(bins = 7, aes(y=..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.7,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of variants arising per transmission") + ylab("Relative frequency") +
  theme(legend.position = "none") 

variants_plot_Treatment <- trans_dists %>% ggplot(aes(x = New_variants)) + 
  geom_histogram(bins = 7, aes(y=2*..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.7,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of variants arising per transmission") + ylab("Relative frequency") +
  theme(legend.position = "none") + facet_wrap(~Treatment)

SNPs_boxplot <- trans_dists %>% ggplot(aes(x = Treatment, y = Dist)) + 
  geom_boxplot(color = "magenta", fill = "blue",alpha= 0.7) +
  ylab("# of SNPs") + xlab("Treatment") +
  theme(legend.position = "none")

variants_boxplot <- trans_dists %>% ggplot(aes(x = Treatment, y = New_variants)) + 
  geom_boxplot(color = "magenta", fill = "blue",alpha= 0.7) +
  ylab("# of New variants") + xlab("Treatment") +
  theme(legend.position = "none")

all_dists_plot <- dists %>% ggplot(aes(x = Dist)) + 
  geom_histogram(bins = 7, aes(y=..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.5,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of SNPs difference") + ylab("Relative frequency") +
  theme(legend.position = "none")


rates_plot <- plot_grid(all_dists_plot,SNPs_plot,variants_plot, nrow = 3, labels = c("A","B","C"))

save_plot("SNP_accumulation.pdf",rates_plot , base_height = 8, base_width = 4)


###Run wilcoxon test to compare
N_SNPs <- trans_dists %>% filter(Treatment == "N") %>% pull(Dist)
W_SNPs <- trans_dists %>% filter(Treatment == "W") %>% pull(Dist)
N_newVar <- trans_dists %>% filter(Treatment == "N") %>% pull(New_variants)
W_newVar <- trans_dists %>% filter(Treatment == "W") %>% pull(New_variants)

wilcox.test(N_newVar,W_newVar)

rates_plot <- plot_grid(all_dists_plot,SNPs_plot,variants_plot, nrow = 3, labels = c("B","C","D"))
SNPs_panel_plot <- plot_grid(Chain_Plots,rates_plot,rel_widths = c(3,1.5), labels = c("A",""))

save_plot("../Images/SNP_acumulation_panel2.pdf", SNPs_panel_plot, base_height = 8, base_width = 11)

### Get average pairwise distance
all_SNPs <- dists %>% pull(SNP_distance)
all_SNPs_fit <-fitdist(all_SNPs,"pois")
print(all_SNPs_fit$estimate)

### Fit the SNPs and New Variants to poisson distribution


SNPs <- trans_dists %>% pull(Dist)
SNPs_fit <-fitdist(SNPs,"pois")
print(SNPs_fit$estimate)
SNPs_lambda <- 0.09424084
######################################
#lambda
#0.09424084 
######################################

Variants <- trans_dists %>% pull(New_variants)
Variants_fit <-fitdist(Variants,"pois")
print(Variants_fit$estimate)
Variants_lambda <- 0.8272251
######################################
#lambda
#0.8272251 
######################################

###Calculate mutation rate with 10000 random draws for time and # mutations
time <- rnorm(10000,mean = 7,sd = 1)
SNPs_reps <- rpois(10000,SNPs_lambda)
Variant_reps <- rpois(10000,Variants_lambda)

#### Mutation rate = SNPS/time*365
SNPs_peryear <- SNPs_reps*365/time
mean(SNPs_peryear)

Variants_peryear <- Variant_reps*365/time
mean(Variants_peryear)

ggplot()+geom_density(aes(Variants_peryear), bw = 25)

ggplot() + geom_density(aes(rnorm(10000,mean = 7,sd = 1)), bw =1,fill = "blue",alpha= 0.8)


##---------------------------------------------------------------------------
# For each comparison count the number of shared hSNPs
##---------------------------------------------------------------------------
for (i in 1: nrow(alleles)){
  if(alleles$Allele_Ratio[i] > 0.025){
    alleles$SNP_type[i] <- "H"
  }
  if(alleles$Allele_Ratio[i] > 0.9){
    alleles$SNP_type[i] <- "Z"
  }
}

for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  hshares <- 0
  weights <- 0
  
  for (qpos in Variant_sites ) {
    q1 <- which(alleles$Pos == qpos & alleles$ID == ID1 )
    q2 <- which(alleles$Pos == qpos & alleles$ID == ID2 )
    qsite <- which(hSites$Pos == qpos)
    freq <- hSites$frequency[qsite]
    weight <- 1/(freq+1)
    if(alleles$Allele_Ratio[q1] > 0.025 & alleles$Allele_Ratio[q2] > 0.025){
      hshares <- hshares + 1
      weight <- weights + weight
    }
  }
  dists$hMatches[i] <- hshares
  dists$hscores[i] <- weights
}



for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  hmatches <- 0
  weighted_matches <- 0
  
  for (qpos in Variants_not_rare) {
  q1 <- which(alleles$Pos == qpos & alleles$ID == ID1 )
  q2 <- which(alleles$Pos == qpos & alleles$ID == ID2 )
  #qsite <- which(hSites$Pos == qpos)
  #freq <- hSites$frequency[qsite]
  #weight <- 1/(freq+1)
  
    if(alleles$SNP_type[q1] %in% c("H","Z") & alleles$SNP_type[q2] %in% c("H","Z")){
    hmatches <- hmatches + 1
    weighted_matches <- weighted_matches + weight
    }
  }
  dists$hMatches[i] <- hmatches
  dists$Weighted_Matches[i] <- weighted_matches
}

##---------------------------------------------------------------------------
# For each timestep, fit hSNPS to poisson and keep average
##---------------------------------------------------------------------------
for (i in 1:40){
  tmp <- trans_dists %>% pull(hMatches)
  
  tmp_fit <-fitdist(tmp,"pois")
  print(tmp_fit$estimate %>% unname())
}


tmp_fit$estimate[1]

## Copy output to bb edit and reformat into a vector
#lambdas <- c(2.973822 ,2.946903 ,2.965909 ,2.782456 ,2.553459 ,2.512968 ,2.459677 ,2.333333 ,2.301663 ,2.264108 ,2.250526 ,2.295146 ,2.263254 ,2.205128 ,2.262136 ,2.165899 ,2.153959 ,2.172222 ,2.159946 ,2.116279 ,2.149686 ,2.191981 ,2.175294 ,2.182156 ,2.165138 ,2.148501 ,2.200577 ,2.205207 ,2.213592 ,2.211409)
cleanhSNP_lambdas <- c(0.3455497,0.4336283,0.3825758,0.3403509,0.2987421,0.2276657,0.1747312,0.1550388,0.1021378,0.07223476,0.05052632,0.03495146,0.02010969,0.02051282,0.01294498,0.01228879,0.008797654,0.006944444,0.00672043,0.003875969)

cleanhZ_lambdas <- c(1.120419,1.030973,0.8787879,0.7719298,0.6477987,0.5244957,0.4381721,0.3979328,0.2992874,0.248307,0.2063158,0.1728155,0.1407678,0.1299145,0.1100324,0.1013825,0.08944282,0.08472222,0.07795699,0.06976744)

plot(x = 1:20, y = log(cleanhSNP_lambdas), xlab="# of transmission steps", ylab="Mean number of shared variants (lambda)", main="Shared variants and transmission", pch=20, col="blue")

Shared_hSNPs_lambdas_plot <- ggplot() + 
  geom_smooth(aes(x = 1:20, y = cleanhSNP_lambdas),se = TRUE , formula= y~x) + 
  geom_point(aes(x = 1:20, y = cleanhSNP_lambdas),pch=20, col="blue") +
  xlab("# of transmission steps") + ylab("Mean shared variants (lambda)") + theme_classic()

Shared_hZ_lambdas_plot <- ggplot() + 
  geom_smooth(aes(x = 1:20, y = cleanhZ_lambdas),se = TRUE , formula= y~x) + 
  geom_point(aes(x = 1:20, y = cleanhZ_lambdas),pch=20, col="blue") +
  xlab("# of transmission steps") + ylab("Mean shared variants (lambda)") + theme_classic()

shared_hSNPs_histogram <- dists %>% filter(timestep < 21) %>% ggplot (aes(Clean_hMatches)) + 
  geom_histogram(bins = 3 , alpha = 0.5, aes(fill = timestep, color = timestep)) +
  xlab("# of shared minority variants") + ylab("Frequency") +
  facet_wrap(~ timestep) + theme_bw() + theme(legend.position = "none")

shared_hZ_histogram <- dists %>% filter(timestep < 21) %>% ggplot (aes(hZ_matches)) + 
  geom_histogram(bins = 5 , alpha = 0.5, aes(fill = timestep, color = timestep)) +
  xlab("# of shared minority variants") + ylab("Frequency") +
  facet_wrap(~ timestep) + theme_bw() + theme(legend.position = "none")

lambdas_panel <- plot_grid(Shared_hSNPs_lambdas_plot,Shared_hZ_lambdas_plot,nrow = 2, labels = c("E","F"))

shared_variants_panel <- plot_grid(shared_hZ_histogram,lambdas_panel, ncol = 2, rel_widths = c(2.5,1), labels = c("D",""))

Variants_panel_plot <- plot_grid(predictive_sites_panel,shared_variants_panel, nrow = 2, rel_heights = c(1,2))

save_plot("../Images/Shared_variants_panel3.pdf", Variants_panel_plot, base_height = 11, base_width = 8)
save_plot("../Images/Shared_variants_panel3.png", Variants_panel_plot, base_height = 11, base_width = 8)

dists %>% filter(timestep < 21) %>% ggplot (aes(x = CleanAr_sum)) + 
  geom_density(bw = 0.05 , alpha = 0.5, fill = timestep, color = timestep) +
  xlab("# of shared minority variants") + ylab("Frequency") +
   theme(legend.position = "none")

##---------------------------------------------------------------------------
# For each comparison count the Allelic ratio differences
##---------------------------------------------------------------------------
for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  AR_mismatches <- 0
  AR_diffs <- vector(mode="numeric", length=0)
  CAR_diffs <- vector(mode="numeric", length=0)
  for (qpos in Variants_not_rare) {
    q1 <- which(alleles$Pos == qpos & alleles$ID == ID1 )
    q2 <- which(alleles$Pos == qpos & alleles$ID == ID2 )
    
    if(alleles$Allele_Ratio[q1] != alleles$Allele_Ratio[q2]){
      AR_mismatches <- AR_mismatches + 1
      diff <- alleles$Allele_Ratio[q1] - alleles$Allele_Ratio[q2]
      AR_diffs <- c(AR_diffs,abs(diff))
    }
  }
  dists$AR_mismatches[i] <- AR_mismatches
  dists$AR_mean[i] <- sum(AR_diffs)/AR_mismatches
  dists$Ar_sum[i] <- sum(AR_diffs)
}  

for (i in 1:nrow(dists)){
  if(dists$CleanAR_mismatches[i] == 0){
    dists$CleanAR_mean[i] <- 0
  }
}


for (i in 1:40){
  tmp <- dists %>% filter(timestep == i) %>% pull(CleanAr_sum)
  
  #tmp_fit <-fitdist(tmp,"pois")
  #print(tmp_fit$estimate %>% unname())
  test <- c(test,mean(tmp))
}
test

plot(x = seq(1,40), y = test)

##---------------------------------------------------------------------------------------
# Plot the changes in hSNPs shared and AR varying sites in the context of timestep
##---------------------------------------------------------------------------------------

dists %>% filter(timestep < 31) %>% ggplot (aes(hZ_matches)) + 
  geom_density(bw = 0.8, alpha = 0.5, aes(group = timestep, color = timestep)) +
  xlab("Number of shared variants") + ylab("Frequency density") +
  scale_color_gradient2(low="navy", mid="white", high="red", midpoint = 15)



dists %>% filter(timestep < 31) %>% ggplot (aes(CleanAr_sum)) + 
  geom_density(bw = 0.1, alpha = 0.5, aes(group = timestep, color = timestep)) +
  xlab("Total change in abundance of alternate allele") + ylab("Frequency density") + 
  scale_color_gradient2(low="navy", mid="white", high="red", midpoint = 15)


##---------------------------------------------------------------------------
# Fit the mean AR change of transmission to a truncated normal distribution
##---------------------------------------------------------------------------

SV_trans <- dists %>% filter(timestep == 1) %>% pull(hZ_matches)
SVT_estimate <- fitdist(SV_trans,"pois")
SVT_estimate$estimate %>% unname()

SVT_labmda <- 1.120419

SV_all <- dists  %>% pull(hZ_matches)
SV_all_fit <- fitdist(SV_all,"pois")
SV_all_lambda <- 0.1537542

AR_mean_trans <- dists %>% filter(timestep == 1) %>% pull(CleanAR_mean)

AR_trans_estimate <- fitdist(AR_mean_trans, "norm")
AR_trans_estimate

AR_trans_mean <- 0.1026025
AR_trans_sd <- 0.1346211

AR_mean_all <- dists %>% pull(CleanAR_mean)

AR_all_estimate <- fitdist(AR_mean_all, "norm")
AR_all_estimate

AR_all_mean <- 0.5172414
AR_all_sd <- 0.3077195

ptrans <- mean(dists$timestep == 1)

#### P(trans|parameter) = P(parameter|trans) * P(trans) / P(parameter)

for(i in 1:nrow(dists)){
  ## Step 1: calculate P(trans|#sharedhz sites)
  P_share <- dpois(dists$hZ_matches[i], lambda = SV_all_lambda)
  P_share_trans <- dpois(dists$hZ_matches[i], lambda = SVT_labmda)
  dists$P_trans_share[i] <- P_share_trans * ptrans / P_share
  
  ## Step 2: calculate P(trans|mean change in allelic ratio)
  P_ARC <- dnorm(dists$CleanAR_mean[i], mean = AR_all_mean, sd = AR_all_sd)
  P_ARC_trans <- dnorm(dists$CleanAR_mean[i], mean = AR_trans_mean, sd = AR_trans_sd)
  dists$P_trans_ARmean[i] <- P_ARC_trans * ptrans / P_ARC
  
  ## Step 3: calculate P(trans|difference in time between )
}
 dists$trans_AR_SV <- dists$P_trans_ARmean * dists$P_trans_share

dists %>% ggplot(aes(x=as.factor(timestep), y = log(trans_AR_SV))) + geom_boxplot() 

samples <- dists %>% dplyr::select(ID2,Ch2,M2) %>% unique()

### For each sample simulate time since outbreak start
for(i in 1:nrow(samples)){
  samples$days[i] <- mean(replicate(1000, (rnorm(samples$M2[i],mean = 7,sd = 2) %>% sum())))
  }


##### Incorporate time into the pairwise comparisons

for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  qt1 <- which(samples$ID2 == ID1)
  qt2 <- which(samples$ID2 == ID2)
  if(dists$ID1[i] == "ICC180"){
    t1 <- 0
    t2 <- samples$days[qt2]
    dists$time_interval[i] <- abs(t1-t2)
  }else{
    t1 <- samples$days[qt1]
    t2 <- samples$days[qt2]
    dists$time_interval[i] <- abs(t1-t2) 
  }
  dists$P_trans_time[i] <- dnorm(dists$time_interval[i], mean = 7, sd = 3)
}

samples %>% filter(M2 == 12) %>% ggplot(aes(x = days)) + geom_density()
dists %>% ggplot(aes(x = as.factor(time_interval), y = log (trans_AR_SV * P_trans_time))) + geom_boxplot()


plot_Ptrans <- function(a){
    dists %>% filter(ID1 ==a | ID2 == a) %>% ggplot(aes(x=time_interval, y = trans_AR_SV, fill = timestep)) + 
    geom_point(aes(color = timestep)) + ylab("P(Transmission)") + xlab("Time interval (days)") + ggtitle(a)
}

plot_Ptrans("N5P12")

Plot1 <- plot_Ptrans("W1P02")

##### For each sample need to rank the hits by 

get_top_hits <- function(a){
  query <- a
  tmp <- dists %>% filter(ID1 %in% query | ID2 %in%query)
}

##---------------------------------------------------------------------------
# Estimate bottleneck sizes using simple normal distribution method
##---------------------------------------------------------------------------


BottleneckLikelihood <- function(c,a,b){
  donor <- a
  recep <- b
  Log_estimate_sum <- 0
  NT <- c
  for (i in c(cleanhSNPS,Zsnps)){
    d <- which((alleles$Pos == i) & (alleles$ID == donor))
    r <- which((alleles$Pos == i) & (alleles$ID == recep))
    
    ARD <- alleles$Allele_Ratio[d]
    ARR <- alleles$Allele_Ratio[r]
    
    if (ARD > 0.01 & ARD < 0.99){
      Variance <- ARD*(1-ARD)/NT
      Log_estimate_sum <- Log_estimate_sum + dnorm(ARR, mean = ARD, sd = Variance, log = TRUE)
    }
  }
  return(-Log_estimate_sum)
}

BottleneckLikelihood(2,"N4P03","N4P04")
trans_dists <- dists  %>% filter(timestep == 1)

for (i in 1:nrow(trans_dists)){
  ID1 <- levels(trans_dists$ID1[i])[trans_dists$ID1[i]]
  ID2 <- levels(trans_dists$ID2[i])[trans_dists$ID2[i]]
  est <- optim(par=1,  fn = BottleneckLikelihood,
                 a = ID1, b = ID2, lower = 1, upper = 100, method = "Brent")
  trans_dists$Nb_Estimate[i] <- as.integer(est$par)
  trans_dists$Nb_loglik[i] <- est$value
}

Nbplot <- trans_dists %>% ggplot(aes(x = M1, y = Nb_Estimate, fill = Ch2, color = Ch2)) + geom_point() + 
          geom_line() + ylim(0,10) + facet_wrap(~Ch2)

plot_grid(Nbplot,Chain_Plots)

a <- "N4P03"
b <- "N4P04"
c <- 2

##---------------------------------------------------------------------------
# Estimate bottleneck sizes using stochastic approach
##---------------------------------------------------------------------------

Stochastic_bottleneck <- function(c,a,b){
  donor <- a
  recep <- b
  Log_estimate_sum <- 0
  NT <- as.integer(c)
  range <- NT-1
  for (i in c(cleanhSNPS,Zsnps)){
    d <- which((alleles$Pos == i) & (alleles$ID == donor))
    r <- which((alleles$Pos == i) & (alleles$ID == recep))
    
    ARD <- alleles$Allele_Ratio[d]
    ARR <- alleles$Allele_Ratio[r]
    if(0.01 < ARD & ARD < 0.99){
      for (k in 1:range){
        binom <- pbinom(k, size=NT, prob=ARD)
        betap <- pbeta(ARR, k, NT-k, ncp = 0, lower.tail = TRUE, log.p = FALSE)
        product <- binom * betap
        
        Log_estimate_sum <- Log_estimate_sum + product
        
        #print(c(i,NT,k,binom,betap,-Log_estimate_sum))
      }
    }
  }
  return(-Log_estimate_sum)
}

Stochastic_bottleneck(2,"N4P03","N4P04")
ID1 <- levels(trans_dists$ID1[i])[trans_dists$ID1[i]]
ID2 <- levels(trans_dists$ID2[i])[trans_dists$ID2[i]]
for (i in 1:nrow(trans_dists)){
  ID1 <- levels(trans_dists$ID1[i])[trans_dists$ID1[i]]
  ID2 <- levels(trans_dists$ID2[i])[trans_dists$ID2[i]]
  est <- optim(par=1,  fn = Stochastic_bottleneck,
               a = ID1, b = ID2, lower = 1, upper = 100, method = "Brent")
  trans_dists$NbStoch_Estimate[i] <- as.integer(est$par)
  trans_dists$NbStoch_loglik[i] <- est$value
}

##---------------------------------------------------------------------------
# Generate allele frequency files to use Katia's BB estimator R script
##---------------------------------------------------------------------------

trans_dists_ARchanges <- trans_dists %>% filter(AR_mismatches > 0)

for (i in 1 : nrow(trans_dists_ARchanges)){
  ID1 <- levels(trans_dists$ID1[i])[trans_dists$ID1[i]]
  ID2 <- levels(trans_dists$ID2[i])[trans_dists$ID2[i]]
  tmp <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(tmp) <- c("D", "R")
  for (j in c(cleanhSNPS,Zsnps)){
    
    d <- which((alleles$Pos == j) & (alleles$ID == a))
    r <- which((alleles$Pos == j) & (alleles$ID == b))
    
    D<- alleles$Allele_Ratio[d]
    R<- alleles$Allele_Ratio[r]
    tmp <- rbind(tmp, c(D,R))
  }
}
 tmp
 write.table(tmp, "../Bottleneck/test.txt", append = FALSE, sep = "\t",
             row.names = FALSE, col.names = FALSE)
####################################################################################
####################################################################################
####################################################################################
####################################################################################
#################### OLD CODE --- OLD CODE --- OLD CODE ############################
####################################################################################
####################################################################################
####################################################################################
####################################################################################




##---------------------------------------------------------------------------
# Analysis of how the number of shared rare hSNPs changes with timestep
##---------------------------------------------------------------------------
tbltest <- table(dists$timestep, dists$Clean_hMatches) %>% prop.table(1)
write.csv(tbltest,"temp_tbl.csv")
tbltemp  <- read.csv("temp_tbl.csv")
colnames(tbltemp) <- c("time","Share0","Share1","Share2")
long_tbl <- melt(tbltemp, id.vars=c("time"))
colnames(long_tbl) <-  c("time", "Shared variants","value")
long_tbl$`Shared variants` <- gsub("Share","",long_tbl$`Shared variants`)

long_tbl %>% filter(time < 31) %>% 
  ggplot(aes(x = time, y = value, fill = `Shared variants`, colour = `Shared variants`)) + geom_smooth() + 
  xlab("Number of transmission Steps") + ylab("Proportion of comparisons") + 
  scale_fill_manual(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) + 
  theme(legend.position = "top")

disttbl <- dists %>% group_by(Same_chain, Clean_hMatches) %>% summarize(count=n())

disttbl %>% ggplot(aes(x = Same_chain, y = count)) +
 geom_bar(aes(fill=as.factor(Clean_hMatches)),position="fill",  stat="identity") +
  scale_fill_manual(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_fill_discrete(name = "Shared Variants")

disttbl

table(dists$Same_chain)


##---------------------------------------------------------------------------
# For each timestep, fit AR mismatches to poisson and keep average
##---------------------------------------------------------------------------
for (i in 1:30){
  tmp <- dists %>% filter(timestep == i) %>% dplyr::select(hMatches)
  
  tmp_fit <-fitdist(tmp$hMatches,"pois")
  print(tmp_fit$estimate)
}

##---------------------------------------------------------------------------
# Create functions to count pairwise matched character states & changes in AR
##---------------------------------------------------------------------------
dists$CleanAR_meanad <- dists$CleanAR_mean 
for (i in 1: nrow(dists)){
  if(dists$CleanAR_mismatches[i] == 0){
    dists$CleanAR_meanad[i] <- 0.0001
  }
}
min(dists$CleanAR_meanad)

##---------------------------------------------------------------------------
# Identify the distribution of AR_mean and fit to distribution
##---------------------------------------------------------------------------
plotdist( dists$CleanAR_meanad, histo = TRUE, demp = TRUE)
descdist( dists$CleanAR_meanad, boot = 1000)

tmp <- AR_change_tbl$AR_diff_mean_adj
ARSum_fit_gamma<- fitdist(tmp,"beta")
summary(ARSum_fit_gamma)

##---------------------------------------------------------------------------
# Create functions to count pairwise matched character states & changes in AR
##---------------------------------------------------------------------------

CS_goodfit <- goodfit(dists$hMatches[which(dists$timestep==7)], "poisson")
Cs_rg <- rootogram(CS_goodfit, xlab = "Number of character state changes")

Get_matches <- function(a,b,c){
  #Initiate variables: count of matches, query indeces (qi)
  nmatches <- 0
  qi1 <- ""
  qi2 <- ""
  qid1 <- levels(a)[a]
  qid2 <- levels(b)[b]
  for (qpos in c) {
    qi1 <- which(AL$Pos == qpos & AL$ID == qid1 )
    qi2 <- which(AL$Pos == qpos & AL$ID == qid2 )
    
    if (AL$SNP_type[qi1] == AL$SNP_type[qi2]) {
      nmatches <- nmatches + 1
    }
  }
  nmatches
}


Get_mismatches <- function(a,b,c){
  #Initiate variables: count of matches, qmuery indeces (qmi)
  nmismatches <- 0
  qmi1 <- ""
  qmi2 <- ""
  qmid1 <- levels(a)[a]
  qmid2 <- levels(b)[b]
  for (qmpos in c) {
    qmi1 <- which(AL$Pos == qmpos & AL$ID == qmid1 )
    qmi2 <- which(AL$Pos == qmpos & AL$ID == qmid2 )
    
    if (AL$SNP_type[qmi1] != AL$SNP_type[qmi2]) {
      nmismatches <- nmismatches + 1
    }
  }
  nmismatches
}

Count_AR_changes <- function(a,b,c){
  #Initiate variables: count of matches, qmuery indeces (qmi)
  nmismatches <- 0
  qmi1 <- ""
  qmi2 <- ""
  qmid1 <- levels(a)[a]
  qmid2 <- levels(b)[b]
  for (qmpos in c) {
    qmi1 <- which(AL$Pos == qmpos & AL$ID == qmid1 )
    qmi2 <- which(AL$Pos == qmpos & AL$ID == qmid2 )
    
    if (AL$Allele_Ratio[qmi1] != AL$Allele_Ratio[qmi2]) {
      nmismatches <- nmismatches + 1
    }
  }
  nmismatches
}

Get_ARdiff<- function(a,b,c){
  #Initiate variables: count of matches, qruery indeces (qri)
  AR_diff <- numeric(0)
  qri1 <- ""
  qri2 <- ""
  qrid1 <- levels(a)[a]
  qrid2 <- levels(b)[b]
  for (qrpos in c) {
    qri1 <- which(AL$Pos == qrpos & AL$ID == qrid1 )
    qri2 <- which(AL$Pos == qrpos & AL$ID == qrid2 )
    tmp <- AL$Allele_Ratio[qri1] - AL$Allele_Ratio[qri2]
    
    AR_diff <- c(AR_diff,abs(tmp))
  }
  sum(AR_diff)
}

  Get_mismatches("ICC180","W1P19",sites)
##---------------------------------------------------------------------------
# Use functions above to update SNP-dists table.
##---------------------------------------------------------------------------
for (i in 1:nrow(dists)){
  dists$Matched_sites[i] <- Get_matches(dists$ID1[i],dists$ID2[i])
  dists$Mismatched_sites[i] <- Get_mismatches(dists$ID1[i],dists$ID2[i])
  dists$Sum_AR_change[i] <- Get_ARdiff(dists$ID1[i],dists$ID2[i])
  dists$AR_change[i] <- Count_AR_changes(dists$ID1[i],dists$ID2[i])
  
}


##---------------------------------------------------------------------------
# Add column for matched chains
##---------------------------------------------------------------------------
for (i in 1:nrow(dists)){
  ch1 <- levels(dists$Ch1[i])[dists$Ch1[i]]
  ch2 <- levels(dists$Ch2[i])[dists$Ch2[i]]
  if (ch1 == ch2){
    dists$Same_chain[i] <- "Yes"
    
  }else {
    dists$Same_chain[i] <- "No"
  }
}

##---------------------------------------------------------------------------
# Calculate the transmission time steps between each pair
##---------------------------------------------------------------------------
for (i in 1:nrow(dists)){
  ch1 <- levels(dists$Ch1[i])[dists$Ch1[i]]
  ch2 <- levels(dists$Ch2[i])[dists$Ch2[i]]
  if (ch1 == ch2){
    dists$timestep[i] <- abs(dists$M1[i] - dists$M2[i])
    
  }else {
    dists$timestep[i] <- dists$M1[i] + dists$M2[i]
  }
}

#Calculate the average AR_change as the SUM_AR/#sites changing AR
for (i in 1:nrow(dists)){
  if (dists$AR_change[i] == 0){
    dists$AR_diff_mean[i] <- 0
  }else{
    dists$AR_diff_mean[i] <- dists$Sum_AR_change[i] / as.numeric(dists$AR_change[i])
  }
}  




##---------------------------------------------------------------------------
#Use set of three complete chains to get estimates on AR_Diff 
##---------------------------------------------------------------------------

Complete_chains <- c("ICC180","N1","N4","W1")
Comp_AL <- AL %>% filter(Chain %in% Complete_chains)

for (i in 1:nrow(Comp_AL)){
  qPos <- Comp_AL$Pos[i]
  qchainlevel <- Comp_AL$Chain[i] 
  qchain <- levels(Comp_AL$Chain[i])[qchainlevel]
  qmouse <- Comp_AL$Mouse_No[i]
  qtypelevel <- Comp_AL$SNP_type[i]
  qtype <- levels(Comp_AL$SNP_type[i])[qtypelevel]
  
  if (qmouse > 1){
    j <- which((Comp_AL$Pos == qPos) & (Comp_AL$Chain == qchainlevel) & (Comp_AL$Mouse_No == qmouse - 1))
    k <- which((Comp_AL$Pos == qPos) & (Comp_AL$Chain == qchain) & (Comp_AL$Mouse_No < qmouse))
    Comp_AL$Priortype[i] <- levels(Comp_AL$SNP_type[j])[Comp_AL$SNP_type[j]]
    Comp_AL$PriorAR[i] <- Comp_AL$Allele_Ratio[j]
    Comp_AL$Type_in_chain[i] <- Comp_AL$SNP_type[i] %in% Comp_AL$SNP_type[k]
  }
  else if (qmouse == 1) {
    j <- which((Comp_AL$Pos == qPos) & (Comp_AL$Chain == "ICC180"))
    Comp_AL$Priortype[i] <- levels(Comp_AL$SNP_type[j])[Comp_AL$SNP_type[j]]
    Comp_AL$PriorAR[i] <- Comp_AL$Allele_Ratio[j]
    Comp_AL$Type_in_chain[i] <- Comp_AL$SNP_type[i] %in% Comp_AL$SNP_type[j]
  }
}

Comp_AL$ARDiff <- Comp_AL$Allele_Ratio - Comp_AL$PriorAR

Comp_AL %>% ggplot(aes(abs(ARDiff))) + 
  geom_density(bw = 0.01, alpha = 0.5) +
  facet_wrap(~ Pos)

##---------------------------------------------------------------------------
#From the training set look at #mismatches at each step
##---------------------------------------------------------------------------

Comp_AL$State_match <- Comp_AL$SNP_type == Comp_AL$Priortype
State_change_tbl <- table(Comp_AL$ID,Comp_AL$State_match)
State_change_tbl <- as.data.frame(State_change_tbl)
State_change_tbl <- spread(State_change_tbl, Var2, Freq)
State_change_tbl <- `colnames<-`(State_change_tbl, c("ID","Mismatched","Matched"))
State_change_tbl$Total <- State_change_tbl$Mismatched + State_change_tbl$Matched 

State_change_tbl <- State_change_tbl %>% filter(!Total == 0)
State_change_tbl <- State_change_tbl %>% filter(!ID == "ICC180")

#Tabulate how many character state changes occur in a given step and plot
table(State_change_tbl$Mismatched)
SC_plot <- qplot(State_change_tbl$Mismatched, geom="histogram", 
      bins = 4, fill=I("blue"), 
      col=I("red"),alpha=I(.2)) +
  xlab("Number of character state changes") + ylab("Frequency") +
  theme_minimal()

#Tabulate the frequency changes in the Allele ratio and plot frquency distribution of AR_Diff
Comp_AL$AR_change <- Comp_AL$ARDiff != 0
AR_change_tbl <- table(Comp_AL$ID,Comp_AL$AR_change)
AR_change_tbl <- as.data.frame(AR_change_tbl)
AR_change_tbl <- spread(AR_change_tbl, Var2, Freq)
AR_change_tbl <- `colnames<-`(AR_change_tbl, c("ID","None","Change"))
AR_change_tbl$Total <- AR_change_tbl$None + AR_change_tbl$Change 

AR_change_tbl <- AR_change_tbl %>% filter(!Total == 0)
AR_change_tbl <- AR_change_tbl %>% filter(!ID == "ICC180")

#Plot frequency changes in sites of AR_Diff for one transmission step
ARC_plot <- qplot(AR_change_tbl$Change, geom="histogram", 
                 bins = 4, fill=I("blue"), 
                 col=I("red"),alpha=I(.2)) +
  xlab("Number of sites changing allele ratio") + ylab("Frequency") +
  theme_minimal()

##---------------------------------------------------------------------------
## Calculate the sum of AR_difference per transmission step
##---------------------------------------------------------------------------

for (i in 1:nrow(AR_change_tbl)){
  qid <- AR_change_tbl$ID[i]
 tmp <- filter(Comp_AL, ID == qid)
 AR_change_tbl$Sum_AR[i] <- sum(abs(tmp$ARDiff))
}

# Calculate the per variable site AR change 

for (i in 1:nrow(AR_change_tbl)){
  if (AR_change_tbl$Change[i] == 0){
    AR_change_tbl$AR_diff_mean[i] <- 0
  }else{
    AR_change_tbl$AR_diff_mean[i] <- AR_change_tbl$Sum_AR[i] / as.numeric(AR_change_tbl$Change[i])
  }
}
# Plot magnitude of AR_Diff per site
AR_sum_plot <- AR_change_tbl %>% ggplot(aes(x = Sum_AR)) + 
   geom_density(bw = 0.05, alpha = 0.2,fill=("blue"), col=("red")) +
  xlab("Magnitude of total allele ratio changes per step") +
  ylab("Probability density") + theme_minimal() 

AR_sum_hist <- qplot(abs(AR_change_tbl$AR_diff_mean_adj), geom="histogram", 
                     binwidth = 0.05,fill=I("blue"),alpha=I(.5)) + 
  xlab("Allele ratio changes per variable site") +
  ylab("Frequency") + theme_minimal()

#Plot values for AR_change per variable site
qplot(abs(AR_change_tbl$Sum_AR/AR_change_tbl$Change), geom="histogram", 
                     binwidth = 0.05,fill=I("blue"),alpha=I(.5)) + 
  xlab("Magnitude of total allele ratio changes per step") +
  ylab("Frequency") + theme_minimal()

AR_change_tbl %>% ggplot(aes(x = Sum_AR/Change)) + 
  geom_density(bw = 0.25, alpha = 0.2,fill=("blue"), col=("red")) +
  xlab("Magnitude of total allele ratio changes per step") +
  ylab("Probability density") + theme_minimal()

# Plot a grid of three graphs
plot_grid(SC_plot, ARC_plot, AR_sum_hist, AR_sum_plot, labels = "AUTO") 


##---------------------------------------------------------------------------
# Fit a poisson distribution to changes is character state and AR.
##---------------------------------------------------------------------------

CS_goodfit <- goodfit(State_change_tbl$Mismatched, "poisson")
Cs_rg <- rootogram(CS_goodfit, xlab = "Number of character state changes")

CS_fit_pois<-fitdist(State_change_tbl$Mismatched,"pois")
coef(CS_fit_pois)

AR_goodfit <- goodfit(AR_change_tbl$Change,"poisson")
AR_rg <- rootogram(AR_goodfit, xlab = "Number of sites changing allele ratio")

AR_fit_pois<-fitdist(AR_change_tbl$Change,"pois")
coef(AR_fit_pois)

plot_grid(Cs_rg,AR_rg, labels = "AUTO")

#Remove values that are 0 or 1 by adding or subtracting 0.00001
for (i in 1:nrow(AR_change_tbl)){
  if(AR_change_tbl$AR_diff_mean[i] == 0){
    AR_change_tbl$AR_diff_mean_adj[i] <- AR_change_tbl$AR_diff_mean[i] + 0.00001
  }
  else if(AR_change_tbl$AR_diff_mean[i] == 1){
    AR_change_tbl$AR_diff_mean_adj[i] <- AR_change_tbl$AR_diff_mean[i] - 0.00001
  }
  else{AR_change_tbl$AR_diff_mean_adj[i] <- AR_change_tbl$AR_diff_mean[i]}
}

##---------------------------------------------------------------------------
# Identify the distribution of SR_sum and fit to distribution
##---------------------------------------------------------------------------
plotdist(AR_change_tbl$AR_diff_mean_adj, histo = TRUE, demp = TRUE)
descdist(AR_change_tbl$AR_diff_mean_adj, boot = 1000)
##---------------------------------------------------------------------------
## Results summary of descdist
## min:  0   max:  1.11 
## median:  0.1675 
## mean:  0.2130455 
## estimated sd:  0.2298165 
## estimated skewness:  1.596145 
## estimated kurtosis:  6.345427 
##---------------------------------------------------------------------------

plotdist(AR_change_tbl$AR_diff_mean, histo = TRUE, demp = TRUE)
descdist(AR_change_tbl$AR_diff_mean, boot = 1000)



tmp <- AR_change_tbl$AR_diff_mean_adj
ARSum_fit_gamma<- fitdist(tmp,"beta")
summary(ARSum_fit_gamma)
##################################
###        Estimate Std. Error
### Shape 1 0.339     0.047
### Shape 2 2.884     0.645
##################################

dbeta(0,0.3391346,2.8836642)

##----------------------------------------------------------------------------------------------------------
# Use the distributions above to calculate probability of observing transmission for full data
##----------------------------------------------------------------------------------------------------------

## P(Charater_changes| Pois(lambda = 0.6212121))
lambdaCS <- 0.6212121
for (i in 1:nrow(dists)){
  k <- dists$Mismatched_sites[i]
  dists$P_CS[i] <- dpois(k,lambdaCS)
}


## P(AR_changes| Pois(lambda = 1.757576 ))
lambdaAR <- 1.757576 
for (i in 1:nrow(dists)){
  k <- dists$AR_change[i]
  dists$P_ARC[i] <- dpois(k,lambdaAR)
}

## P(Per_site_ARchange | Beta(0.3391346,2.8836642 ))
for (i in 1:nrow(dists)){
  if(dists$AR_diff_mean[i] == 0){
    dists$AR_diff_mean_adj[i] <- dists$AR_diff_mean[i] + 0.00001
  }
  else if(dists$AR_diff_mean[i] == 1){
    dists$AR_diff_mean_adj[i] <- dists$AR_diff_mean[i] - 0.00001
  }
  else{dists$AR_diff_mean_adj[i] <- dists$AR_diff_mean[i]}
}

for (i in 1:nrow(dists)){
  dists$P_ARDiff[i] <- dbeta(dists$AR_diff_mean_adj[i],0.3391346,2.8836642)
}

## Calculate the composite probability as the product of probabiloity of observing character states and AR difference per variable site

dists$P_trans <- dists$P_CS * dists$P_ARDiff


scores<- numeric(0)
for(k in 0:8){
score <- (lambdaCS^k)*(exp(lambdaCS))/factorial(k)
scores <- c(scores,score)
}
dpois(0, 0.6212121)

##----------------------------------------------------------------------------------------------------------
# Calculate P_sig and P-Pi from full data
##----------------------------------------------------------------------------------------------------------
Pi_goodfit <- goodfit(dists$Mismatched_sites, "poisson")
Pi_rg <- rootogram(Pi_goodfit, xlab = "Number of character state changes")

Pi_fit_pois<-fitdist(dists$Mismatched_sites,"pois")
coef(Pi_fit_pois)


############
#lambda 
#2.561645 
##########

lambdaPr <- 2.561645
for (i in 1:nrow(dists)){
  k <- dists$Mismatched_sites[i]
  dists$P_Pi[i] <- dpois(k,lambdaPr)
}

## 
plotdist(dists$AR_diff_mean_adj, histo = TRUE, demp = TRUE)
descdist(dists$AR_diff_mean_adj, boot = 1000)
# Uniform distribution
dunif(0.3, min = 0, max = 1, log = FALSE)
runif(20901, min = 0, max = 1)
dunif(dists$AR_diff_mean_adj, a = 0.0001, b = 1)
##----------------------------------------------------------------------------------------------------------
# Use probability density to calculate plot heatmaps
##----------------------------------------------------------------------------------------------------------

ggplot(data = dists) + geom_tile(aes (x = ID2, y = ID1 , fill = Sum_AR_change))
ggplot(data = dists) + geom_tile(aes (x = ID1, y = ID2 , value = Dist))

dists %>% filter(Dist==0) %>% 
  ggplot(aes(x=as.factor(timestep), y = log(P_trans))) + 
  geom_violin(fill='#A4A4A4', color="blue") + xlab("# of transmission steps") + ylab("Probability of transmission (log)") +
  theme_minimal()
dists %>% 
  ggplot(aes(x=as.factor(timestep), y = Sum_AR_change)) + 
  geom_boxplot(fill='#A4A4A4', color="blue") + xlab("# of transmission steps") + ylab("Probability of transmission (log)") +
   theme_minimal()

dists %>% 
  ggplot(aes(x=as.factor(timestep), y = Dist)) + 
  geom_violin(fill='#A4A4A4', color="blue") + xlab("# of transmission steps") + ylab("Probability of transmission (log)") +
  theme_minimal()
dists %>% filter(Dist < 2) %>% 
  ggplot(aes(x=as.factor(timestep), y = log(P_trans))) + 
  geom_boxplot() + xlab("# of transmission steps") + ylab("Probability of transmission (log)") +
  theme_minimal()


##----------------------------------------------------------------------------------------------------------
# Test the metric on chains N3, N5 and W4
##----------------------------------------------------------------------------------------------------------

Ch_subset <- c("ICC180","N3","N5","W4")

dist_subset <- dists %>% filter(Ch1 %in% Ch_subset & Ch2 %in% Ch_subset)

ggplot(data = dist_subset) + geom_tile(aes (x = ID2, y = ID1 , fill = log(P_trans)))

tmp_dist <- dists %>% filter((Dist == 0) & (ID1 == "N3P12" | ID2 == "N3P12"))

tmp_dist %>% dplyr::select(ID1, ID2,"timestep", "P_trans")

tmp_dist %>% 
  ggplot(aes(x=as.factor(timestep), y = log(P_trans))) + 
  geom_point() + xlab("# of transmission steps") + ylab("Probability of transmission (log)") +
  main("An example with strain N3P12 and identical genomes") + theme_minimal()


old_dists <- read.csv("../Transmission_likelihood/Pairwise_comparisons.csv")

plot_sample <- function(a){
  tmp_dist <- old_dists %>% filter((ID1 == a | ID2 == a))
  
  #tmp_dist %>% dplyr::select(ID1, ID2,"timestep", "P_trans")
  
  p1 <- tmp_dist %>% 
    ggplot(aes(x=transmission_steps, y = Mean_allelic_ratio_change)) + 
    geom_point(aes(color = Same_chain)) + xlab("# of transmission steps") + 
    ylab("Mean AR Change")  +
     theme_minimal() + theme(legend.position = "none")
  p2 <- tmp_dist %>% 
      ggplot(aes(x=transmission_steps, y = Total_allelic_ratio_change)) + 
    geom_point(aes(color = Same_chain)) + xlab("# of transmission steps") + 
    ylab("Total AR Change") + ggtitle(a) + 
    theme_minimal() + theme(legend.position = "blank")
  p3 <- tmp_dist %>% 
    ggplot(aes(x=transmission_steps, y = P_trans_ARmean)) + 
    geom_point(aes(color = Same_chain)) + xlab("# of transmission steps") + 
    ylab("Likelihood of transmission") + 
    theme_minimal() + theme(legend.position = "none") 
  
  
  
  plot_grid(p1,p2,p3,ncol=3, rel_widths = c(1,1,1))
  
}
plot_sample("N3P12")
plot_sample("W1P12")

##----------------------------------------------------------------------------------------------------------
# Find correlation between SNP Dists, H SNP dists and P_trans with trans_steps
##----------------------------------------------------------------------------------------------------------

cor(dists$P_t, dists$timestep, method = c("pearson", "kendall", "spearman"))

ggscatter(data = (filter (dists, Dist == 0)), x = "timestep", y = "AR_diff_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")

############# Old code ####################################################
#Comp_IDs <- Comp_AL$ID
#k <- ""
#mmcounts <- ""
#for (k in Comp_IDs){
  #print (k)
#  tmp <- sum(Comp_AL$SNP_type != Comp_AL$Priortype & Comp_AL$ID == k)
 # mmcounts <- c(mmcounts,tmp)
#}
############################################################################
##---------------------------------------------------------------------------
#Plot line graphs showing how allele rations change over the course
##---------------------------------------------------------------------------

AR_plot <- AL %>% filter(Pos %in% Zsnps) %>% ggplot(aes(x=Mouse_No, y= Allele_Ratio, fill = Chain, color = Chain )) + geom_line() +
  facet_wrap_paginate(~ Pos) + 
  theme_bw() + xlab("Step in transmission chain") + ylab("Proportion of reads mapped to alternate allele") +
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#fafa07"))

AL %>% filter(Pos %in% Good_sites) %>% ggplot(aes(x=Mouse_No, y= Allele_Ratio, fill = Chain, color = Chain )) + geom_line() + geom_hline(yintercept =0.9) +
  facet_wrap(~ Pos) + xlab("Time (position in transmission Chain)") + ylab("Proportion of reads mapping to alternative allele")

##---------------------------------------------------------------------------
## Plot graphs alongside the phylogenetic tree
##---------------------------------------------------------------------------
ctree <- read.tree("../RAxML_bestTree.test")
data <- read.csv("../Wide_table_Allele_ratios.csv", check.names = FALSE)
p <- ggtree(ctree) %<+% data
p2 <- p + geom_tippoint(aes(shape = Chain, color = Chain, label = Mouse_No)) + 
  scale_shape_manual(values=c(0,1,2,3, 5,6,8,9,16, 17,11)) +
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99"))


plot_grid(p2, AR_plot, labels = "AUTO", rel_widths = c(1, 2))


##---------------------------------------------------------------------------
# Plot a heatmap of the distance matrix
##---------------------------------------------------------------------------

ggplot(data = dists, aes(x=ID2, y=ID1, fill=Dist)) + 
  geom_tile() + theme(axis.text.y = element_text(size=3)) +
  theme(axis.text.x = element_text(size=3, angle = 90))

ggplot(data = dists, aes(x=ID2, y=ID1, fill = Mismatched_sites)) + 
  geom_tile() + theme(axis.text.y = element_text(size=3)) + theme(axis.text.x = element_text(size=3, angle = 90))
boxplot(dists$AR_diff_mean~dists$timestep)
