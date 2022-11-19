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

setwd("~/Projects/C._rodentium/gvcf_reanalysis/Transmission_likelihood/")

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
ICC_hSNPs <- alleles$Pos[which(alleles$ID == "ICC180" & alleles$SNP_type == "H")]
 
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
# For each transmission record the SNPs distance and get average
##---------------------------------------------------------------------------
alleles %>% filter(Pos %in% Zsnps) %>% filter(!ID == "ICC180") %>% 
  ggplot(aes(x=Mouse_No, y= Allele_Ratio, fill = as.factor(Pos), color = as.factor(Pos))) + 
  geom_line() + geom_point(alpha = 0.7) + geom_hline(yintercept =0.9) +
  facet_wrap(~ Chain, ncol = 3) + xlab("Time (position in transmission Chain)") + ylab("Proportion of reads mapping to alternative allele")

N_plots <- alleles %>% filter(Pos %in% c(Zsnps,cleanhSNPS),!ID == "ICC180",Chain %in% c("N1","N2","N3","N4","N5")) %>%
     ggplot(aes(x=Mouse_No, y= Allele_Ratio, fill = as.factor(Pos), color = as.factor(Pos))) + 
     geom_line() + geom_point(alpha = 0.7) + geom_hline(yintercept =0.9) +
     facet_wrap(~ Chain, ncol = 1) + xlab("Time (position in transmission Chain)") + 
     ylab("Proportion of reads mapping to alternative allele") + 
     theme(legend.position = "none")

W_plots <- alleles %>% filter(Pos %in% c(Zsnps,cleanhSNPS),!ID == "ICC180",!Chain %in% c("N1","N2","N3","N4","N5")) %>%
     ggplot(aes(x=Mouse_No, y= Allele_Ratio, fill = as.factor(Pos), color = as.factor(Pos))) + 
     geom_line() + geom_point(alpha = 0.7) + geom_hline(yintercept =0.9) +
     facet_wrap(~ Chain, ncol = 1) + xlab("Time (position in transmission Chain)") + ylab("") +
     theme(axis.text.y = element_blank(), legend.title = element_blank())+ guides(fill=guide_legend(ncol=1))

leg <- as_ggplot(get_legend(W_plots)) 

W_plot <- W_plots + theme(legend.position = "none")

Chain_Plots <- plot_grid(N_plots, W_plots, ncol=2, rel_widths = c(1,1.25))


##---------------------------------------------------------------------------
# For each comparison count how many times a new variant emerges
##---------------------------------------------------------------------------


for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  new_sites <- 0
  
  for (qpos in c(cleanhSNPS,Zsnps)) {
    q1 <- which(alleles$Pos == qpos & alleles$ID == ID1 )
    q2 <- which(alleles$Pos == qpos & alleles$ID == ID2 )
    
    if(alleles$SNP_type[q1] == "R" & alleles$SNP_type[q2] %in% c("H","Z")){
      new_sites <- new_sites + 1
    }
  }
  dists$New_variants[i] <- new_sites
}  

SNPs_plot <- dists %>% filter(timestep == 1) %>% ggplot(aes(x = Dist)) + 
  geom_histogram(bins = 3, aes(y=..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.5,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of SNPs per transmission") + ylab("Relative frequency") +
  theme(legend.position = "none")

variants_plot <- dists %>% filter(timestep == 1) %>% ggplot(aes(x = New_variants)) + 
  geom_histogram(bins = 4, aes(y=..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.5,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of variants arising per transmission") + ylab("Relative frequency") +
  theme(legend.position = "none")

all_dists_plot <- dists %>% ggplot(aes(x = Dist)) + 
  geom_histogram(bins = 7, aes(y=..ncount../sum(..ncount..)),color = "blue", fill = "blue",alpha= 0.4) + 
  geom_density(bw = 0.5,color = "magenta", fill = "magenta",alpha= 0.7 ) + ylim(0,1) +
  xlab("# of SNPs difference") + ylab("Relative frequency") +
  theme(legend.position = "none")

rates_plot <- plot_grid(all_dists_plot,SNPs_plot,variants_plot, nrow = 3, labels = c("B","C","D"))
SNPs_panel_plot <- plot_grid(Chain_Plots,rates_plot,rel_widths = c(3,1.5), labels = c("A",""))

save_plot("../Images/SNP_acumulation_panel2.pdf", SNPs_panel_plot, base_height = 8, base_width = 11)

### Fit the SNPs and New Variants to poisson distribution
SNPs <- dists %>% filter(timestep == 1) %>% pull(Dist)
SNPs_fit <-fitdist(SNPs,"pois")
print(SNPs_fit$estimate)
SNPs_lambda <- 0.09424084
######################################
#lambda
#0.09424084 
######################################

Variants <- dists %>% filter(timestep == 1) %>% pull(New_variants)
Variants_fit <-fitdist(Variants,"pois")
print(Variants_fit$estimate)
Variants_lambda <- 0.2879581
######################################
#lambda
#0.2879581 
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
# Look at HSNP sites and ask the following question:
# 1) How many isolates have it as a hSNPs
# 2) How frequently is it shared between transmission pairs
# 3) how frequently is it shared between non transmission pairs
##---------------------------------------------------------------------------

hSites <- alleles %>% dplyr::select(Pos) %>% unique()

for(i in 1:nrow(hSites)){
  if(hSites$Pos[i] %in% Zsnps){
    hSites$Designation[i] <- "SNP"
  }else if (hSites$Pos[i] %in% Noisy){
    hSites$Designation[i] <- "Noisy"
  }else if (hSites$Pos[i] %in% cleanhSNPS){
    hSites$Designation[i] <- "hSNPs"
  }else{hSites$Designation[i] <- "Rare"}
}

## Quantifies how many times a site is shared across transmission pairs, cluster pairs and chain pairs
for (i in 1:nrow(hSites)){
  hSites$frequency[i] <-  alleles %>% filter(Pos == hSites$Pos[i], SNP_type %in% c("H","Z")) %>% nrow()
  trans_hmatches <- 0
  no_trans_hmatches <- 0
  trans_hvmatches <- 0
  no_trans_hvmatches <- 0
  
  cluster_hmatches <- 0
  no_cluster_hmatches <- 0
  cluster_hvmatches <- 0
  no_cluster_hvmatches <- 0
  
  chain_hmatches <- 0
  no_chain_hmatches <- 0
  chain_hvmatches <- 0
  no_chain_hvmatches <- 0
  
  for(j in 1:nrow(dists)){
    ID1 <- levels(dists$ID1[j])[dists$ID1[j]]
    ID2 <- levels(dists$ID2[j])[dists$ID2[j]]
    index1 <- which(alleles$Pos == hSites$Pos[i] & alleles$ID == ID1 )
    index2 <- which(alleles$Pos == hSites$Pos[i] & alleles$ID == ID2 )
    
    if(dists$timestep[j] == 1){
      if(alleles$SNP_type[index1] == "H" && alleles$SNP_type[index2] =="H"){
        trans_hmatches <- trans_hmatches + 1
      }
      if(alleles$SNP_type[index1] %in% c("H","Z") && alleles$SNP_type[index2] %in% c("H","Z")){
        trans_hvmatches <- trans_hvmatches + 1
      }
    } 
    if(dists$timestep[j] > 1){
      if(alleles$SNP_type[index1] == "H" && alleles$SNP_type[index2] =="H"){
        no_trans_hmatches <- no_trans_hmatches + 1
      }
      if(alleles$SNP_type[index1] %in% c("H","Z") && alleles$SNP_type[index2] %in% c("H","Z")){
        no_trans_hvmatches <- no_trans_hvmatches + 1
      }
    }
    if(dists$timestep[j] < 6){
      if(alleles$SNP_type[index1] == "H" && alleles$SNP_type[index2] =="H"){
        cluster_hmatches <- cluster_hmatches + 1
      }
      if(alleles$SNP_type[index1] %in% c("H","Z") && alleles$SNP_type[index2] %in% c("H","Z")){
        cluster_hvmatches <- cluster_hvmatches + 1
      }
    } 
    if(dists$timestep[j] > 5){
      if(alleles$SNP_type[index1] == "H" && alleles$SNP_type[index2] =="H"){
        no_cluster_hmatches <- no_cluster_hmatches + 1
      }
      if(alleles$SNP_type[index1] %in% c("H","Z") && alleles$SNP_type[index2] %in% c("H","Z")){
        no_cluster_hvmatches <- no_cluster_hvmatches + 1
      }
    }
    if(dists$Same_chain[j] == "Yes"){
      if(alleles$SNP_type[index1] == "H" && alleles$SNP_type[index2] =="H"){
        chain_hmatches <- chain_hmatches + 1
      }
      if(alleles$SNP_type[index1] %in% c("H","Z") && alleles$SNP_type[index2] %in% c("H","Z")){
        chain_hvmatches <- chain_hvmatches + 1
      }
    } 
    if(dists$Same_chain[j] == "No"){
      if(alleles$SNP_type[index1] == "H" && alleles$SNP_type[index2] =="H"){
        no_chain_hmatches <- no_chain_hmatches + 1
      }
      if(alleles$SNP_type[index1] %in% c("H","Z") && alleles$SNP_type[index2] %in% c("H","Z")){
        no_chain_hvmatches <- no_chain_hvmatches + 1
      }
    }
  }
  
  hSites$trans_shared[i] <- trans_hmatches
  hSites$Notrans_shared[i] <- no_trans_hmatches
  hSites$trans_hvshared[i] <- trans_hvmatches
  hSites$Notrans_hvshared[i] <- no_trans_hvmatches
  
  hSites$cluster_shared[i] <- cluster_hmatches
  hSites$Nocluster_shared[i] <- no_cluster_hmatches
  hSites$cluster_hvshared[i] <- cluster_hvmatches
  hSites$Nocluster_hvshared[i] <- no_cluster_hvmatches
  
  hSites$chain_shared[i] <- chain_hmatches
  hSites$Nochain_shared[i] <- no_chain_hmatches	
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
  hSites$PTS[i] <- (hSites$trans_shared[i]/ntrans)*ptrans/((hSites$trans_shared[i] + hSites$Notrans_shared[i])/nrow(dists))
  hSites$PTCl[i] <- (hSites$cluster_shared[i]/ncluster)*pcluster/((hSites$cluster_shared[i] + hSites$Nocluster_shared[i])/nrow(dists))
  hSites$PTCh[i] <- (hSites$chain_shared[i]/nchain)*pchain/((hSites$chain_shared[i] + hSites$Nochain_shared[i])/nrow(dists))
  
  hSites$PTS_hv[i] <- (hSites$trans_hvshared[i]/ntrans)*ptrans/((hSites$trans_hvshared[i] + hSites$Notrans_hvshared[i])/nrow(dists))
  hSites$PTCl_hv[i] <- (hSites$cluster_hvshared[i]/ncluster)*pcluster/((hSites$cluster_hvshared[i] + hSites$Nocluster_hvshared[i])/nrow(dists))
  hSites$PTCh_hv[i] <- (hSites$chain_hvshared[i]/nchain)*pchain/((hSites$chain_hvshared[i] + hSites$Nochain_hvshared[i])/nrow(dists))
}

PTS_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = PTS, color = Designation, )) + 
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  ylab("P(Transmission|Shared variant)") + xlab("") +
  geom_point() + theme(legend.position='none')

PTCl_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = PTCl, color = Designation, )) + 
  ylab("P(Cluster|Shared variant)") + xlab("Proportion with variant (%)") +
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  geom_point() + theme(legend.position='none')

PTCh_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = PTCh, color = Designation, )) + 
  ylab("P(Same chain|Shared variant)") + xlab("Proportion of strains with variant (%)") +
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  geom_point() + theme(legend.position="right") + theme(legend.title = element_blank())


PTS_hv_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = PTS_hv, color = Designation, )) + 
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  ylab("P(Transmission|Shared)") + xlab("") +
  geom_point() + theme(legend.position='none')

PTCl_hv_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = PTCl_hv, color = Designation, )) + 
  ylab("P(Cluster|Shared)") + xlab("Proportion with variant (%)") +
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  geom_point() + theme(legend.position='none')

PTCh_hv_plot <- hSites %>% ggplot(aes(x = frequency/2.05, y = PTCh_hv, color = Designation, )) + 
  ylab("P(Same chain|Shared)") + xlab("") +
  scale_color_brewer(type = "qual", guide = "legend", palette = "Dark2") + 
  geom_point() + theme(legend.position="right") + theme(legend.title = element_blank())


x_blanker <- theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

legend <- as_ggplot(get_legend(PTCh_hv_plot))

predictive_sites_panel <- plot_grid(PTS_hv_plot,PTCl_hv_plot,PTCh_hv_plot + theme(legend.position='none'), ncol =  3, labels = "AUTO")
predictive_sites_panel <- plot_grid(predictive_sites_panel,legend, rel_widths = c(9,1))
##---------------------------------------------------------------------------
# For each comparison count the number of shared hSNPs
##---------------------------------------------------------------------------


for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  hmatches <- 0
  chmatches <- 0
  zmatches <- 0
  
  for (qpos in hSNPs) {
  q1 <- which(alleles$Pos == qpos & alleles$ID == ID1 )
  q2 <- which(alleles$Pos == qpos & alleles$ID == ID2 )
  
    if(alleles$SNP_type[q1] == "H" & alleles$SNP_type[q2] == "H"){
    hmatches <- hmatches + 1
    }
  }
  for (qcpos in cleanhSNPS) {
    qc1 <- which(alleles$Pos == qcpos & alleles$ID == ID1 )
    qc2 <- which(alleles$Pos == qcpos & alleles$ID == ID2 )
    
    if(alleles$SNP_type[qc1] == "H" & alleles$SNP_type[qc2] == "H"){
      chmatches <- chmatches + 1
    }
  }
  for (qzpos in Zsnps) {
    qz1 <- which(alleles$Pos == qzpos & alleles$ID == ID1 )
    qz2 <- which(alleles$Pos == qzpos & alleles$ID == ID2 )
    
    if(alleles$SNP_type[qz1] %in% c("H","Z") & alleles$SNP_type[qz2] %in% c("H","Z")){
      zmatches <- zmatches + 1
    }
  }
  
  dists$hMatches[i] <- hmatches 
  dists$Clean_hMatches[i] <- chmatches
  dists$hZ_matches[i] <- chmatches + zmatches
}


##---------------------------------------------------------------------------
# For each timestep, fit hSNPS to poisson and keep average
##---------------------------------------------------------------------------
for (i in 1:20){
  tmp <- dists %>% filter(timestep == i) %>% pull(Clean_hMatches)
  
  tmp_fit <-fitdist(tmp,"pois")
  print(tmp_fit$estimate %>% unname())
}

for (i in 1:20){
  tmp <- dists %>% filter(timestep == i) %>% pull(hZ_matches)
  
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

##---------------------------------------------------------------------------
# For each comparison count the Allelic ratio differences
##---------------------------------------------------------------------------
for (i in 1:nrow(dists)) {
  ID1 <- levels(dists$ID1[i])[dists$ID1[i]]
  ID2 <- levels(dists$ID2[i])[dists$ID2[i]]
  AR_mismatches <- 0
  CAR_mismatches <- 0
  AR_diffs <- vector(mode="numeric", length=0)
  CAR_diffs <- vector(mode="numeric", length=0)
  for (qpos in c(hSNPs,Zsnps)) {
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
  
  for (cqpos in c(cleanhSNPS,Zsnps)) {
    cq1 <- which(alleles$Pos == cqpos & alleles$ID == ID1 )
    cq2 <- which(alleles$Pos == cqpos & alleles$ID == ID2 )
    
    if(alleles$Allele_Ratio[cq1] != alleles$Allele_Ratio[cq2]){
      CAR_mismatches <- CAR_mismatches + 1
      Cdiff <- alleles$Allele_Ratio[cq1] - alleles$Allele_Ratio[cq2]
      CAR_diffs <- c(CAR_diffs,abs(Cdiff))
    }
  }
  dists$CleanAR_mismatches[i] <- CAR_mismatches
  dists$CleanAR_mean[i] <- sum(CAR_diffs)/CAR_mismatches
  dists$CleanAr_sum[i] <- sum(CAR_diffs)
  #dists$ARDiffs[i] <- AR_diffs
  #dists$CARDiffs[i] <- CAR_diffs
}  

for (i in 1:nrow(dists)){
  if(dists$CleanAR_mismatches[i] == 0){
    dists$CleanAR_mean[i] <- 0
  }
}

##---------------------------------------------------------------------------------------
# Plot the changes in hSNPs shared and AR varying sites in the context of timestep
##---------------------------------------------------------------------------------------

dists %>% filter(timestep < 31) %>% ggplot (aes(hZ_matches)) + 
  geom_density(bw = 0.8, alpha = 0.5, aes(group = timestep, color = timestep)) +
  xlab("Number of shared variants") + ylab("Frequency density") +
  scale_color_gradient2(low="navy", mid="white", high="red", midpoint = 15)



dists %>% filter(timestep < 31) %>% ggplot (aes(CleanAR_mean)) + 
  geom_density(bw = 0.1, alpha = 0.5, aes(group = timestep, color = timestep)) +
  xlab("Mean change in abundance of alternate allele") + ylab("Frequency density") + 
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

samples %>% filter(M2 == 12) %>% ggplot(aes(x = days)) + geom_density()

get_top_hits <- function(a){
  query <- a
  tmp <- dists %>% filter(ID1 %in% query | ID2 %in%query)
}

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
plot_sample <- function(a){
  tmp_dist <- dists %>% filter((Dist < 2) & (ID1 == a | ID2 == a))
  
  #tmp_dist %>% dplyr::select(ID1, ID2,"timestep", "P_trans")
  
  tmp_dist %>% 
    ggplot(aes(x=as.factor(timestep), y = log(P_trans))) + 
    geom_point(aes(color = Same_chain)) + xlab("# of transmission steps") + ylab("Probability density (log)") +
     theme_minimal()
  
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
