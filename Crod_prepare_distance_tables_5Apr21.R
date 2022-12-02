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

for (i in 1:nrow(hSites)){
  hSites$frequency[i] <-  alleles %>% filter(Pos == hSites$Pos[i], Allele_Ratio > 0.025) %>% nrow()
}

Variant_sites <- hSites %>% filter(!Pos %in% Noisy) %>% filter(frequency < 200, frequency != 0) %>% pull(Pos)
Variants_not_rare <- hSites %>% filter(!Pos %in% Noisy) %>% filter(frequency > 2, frequency < 200) %>% pull(Pos)


## Quantifies how many times a site is shared across transmission pairs, cluster pairs and chain pairs
for (i in 1:nrow(hSites)){
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

write.csv(hSites,"Data_tables/Sites_predictiveness.csv", row.names = FALSE)

##---------------------------------------------------------------------------
# For each comparison count number of shared variants
##---------------------------------------------------------------------------

for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  shared_sites <- 0
  
  for (qpos in Variants_not_rare) {
    q1 <- which(alleles$Pos == qpos & alleles$ID == ID1 )
    q2 <- which(alleles$Pos == qpos & alleles$ID == ID2 )
    
    if(alleles$Allele_Ratio[q1] > 0.025 & alleles$Allele_Ratio[q2] > 0.025 ){
      shared_sites <- shared_sites + 1
    }
  }
  dists$Shared_variants[i] <- shared_sites
}  

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
  if(dists$AR_mismatches[i] == 0){
    dists$AR_mean[i] <- 0
  }
}

##---------------------------------------------------------------------------
# Fit the mean AR change of transmission to a truncated normal distribution
##---------------------------------------------------------------------------

SV_trans <- dists %>% filter(timestep == 1) %>% pull(Shared_variants)
SVT_estimate <- fitdist(SV_trans,"pois")
SVT_estimate$estimate %>% unname()

SVT_labmda <- 1.298429

SV_all <- dists  %>% pull(Shared_variants)
SV_all_fit <- fitdist(SV_all,"pois")
SV_all_fit$estimate
SV_all_lambda <- 0.1713056

AR_mean_trans <- dists %>% filter(timestep == 1) %>% pull(AR_mean)

AR_trans_estimate <- fitdist(AR_mean_trans, "norm")
AR_trans_estimate

AR_trans_mean <- 0.1118154
AR_trans_sd <- 0.1197998

AR_mean_all <- dists %>% pull(AR_mean)

AR_all_estimate <- fitdist(AR_mean_all, "norm")
AR_all_estimate

AR_all_mean <- 0.4733320
AR_all_sd <- 0.2808522

ptrans <- mean(dists$timestep == 1)

#### P(trans|parameter) = P(parameter|trans) * P(trans) / P(parameter)

for(i in 1:nrow(dists)){
  ## Step 1: calculate P(trans|#sharedhz sites)
  P_share <- dpois(dists$Shared_variants[i], lambda = SV_all_lambda)
  P_share_trans <- dpois(dists$Shared_variants[i], lambda = SVT_labmda)
  dists$P_trans_share[i] <- P_share_trans * ptrans / P_share
  
  ## Step 2: calculate P(trans|mean change in allelic ratio)
  P_ARC <- dnorm(dists$AR_mean[i], mean = AR_all_mean, sd = AR_all_sd)
  P_ARC_trans <- dnorm(dists$AR_mean[i], mean = AR_trans_mean, sd = AR_trans_sd)
  dists$P_trans_ARmean[i] <- P_ARC_trans * ptrans / P_ARC
  
  ## Step 3: calculate P(trans|difference in time between )
}

write.csv(dists, "Data_tables/Pairwise_comparisons.csv", row.names = FALSE)




#######/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/#######

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










