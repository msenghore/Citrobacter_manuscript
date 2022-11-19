library(ggplot2)
library(ggforce)
library(plyr)
library("dplyr")
library(tidyr)
library("beeswarm")
library("viridis")
library("fitdistrplus")
library("hrbrthemes")
library("vcd")
library("ggtree")
library(cowplot)
library(tidyverse)
library("ggpubr")
library("ActuDistns")
library("reshape2")

setwd("~/Projects/C._rodentium/gvcf_reanalysis/Bottleneck/")


## Read in the SNP distance table and the full table with allele ratios at each site.
dists <- read.csv("Crod_ID-Ordered_Dist_table.csv")
alleles <- read.csv("All_snps_summary_long_clean.csv")

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

### Identify sites of interest based on previous analysis

#Sites deemed to be noise
Noisy <- c(123360,175076,408983,645840,890802,1086564,1417667,2169284,2172181,2194913,2584041,3040456,4623423,5081862)

#Assign query sites: sites with Zsnps or at least 6 hSNPs
#Remove the noisy sites
Good_sites <- c(2438536,3996563,1285412,2835878,793586,1039375,1229417,17805,1749372,1896011,4987380,492511,2372393,175076,2194913,4623423,645840,3040456,2584041,1417667,123360,890802,937620,2172181,408983,3638848,1531388,2169284,5081862,1086564,4298690,1415559,4854857,1032727,1265896,3928225)
sites <- Good_sites[!Good_sites %in% Noisy]

# Check for overlap
Zsnps <- c(2372393 , 492511 , 1896011 , 4987380 , 17805 , 1749372 , 1229417 , 1039375 , 793586 , 1285412 , 2835878 , 3996563 , 2438536)
Zsnps[Zsnps %in% sites]

###################################################################################
### Look at how the variance changes for different values or AR and NT

NT_AR <- read.csv("Variance_NT_Ar_combinations.csv")
AR_NT_long <- melt(NT_AR, id.vars=c("NT"))

AR_NT_long <- AR_NT_long %>% mutate(variable = str_replace(variable, "X", ""))
AR_NT_long <- rename(AR_NT_long, c("Allele_Ratio"="variable", "Variance"="value"))


ggplot(data=AR_NT_long, aes(x=NT, y= -log(Variance))) +
  geom_line() + facet_wrap(~ Allele_Ratio) 

ggplot(data=AR_NT_long, aes(x=NT, y=Variance)) +
  geom_line() + facet_wrap(~ Allele_Ratio) + xlab("Bottleneck size")
###################################################################################

### Create a new dataframe where the allele ratio is adjusted to eliminate 0's and 1's

alleles_ad <- alleles 
for (i in 1:nrow(alleles_ad)){
  if(alleles_ad$Allele_Ratio[i] == 0){
    alleles_ad$Allele_Ratio_adj[i] <- alleles_ad$Allele_Ratio[i] + 0.0001
  }
  else if(alleles_ad$Allele_Ratio[i] == 1){
    alleles_ad$Allele_Ratio_adj[i] <- alleles_ad$Allele_Ratio[i] - 0.0001
  }
  else{alleles_ad$Allele_Ratio_adj[i] <- alleles_ad$Allele_Ratio[i]}
}

### Write a function to take two genomes, calculate likelihood for bottleneck size

BottleneckLikelihood <- function(a,b,c){
  donor <- a
  recep <- b
  Log_estimate_sum <- 0
  NT <- c
  for (i in sites){
      d <- which((alleles$Pos == i) & (alleles$ID == donor))
      r <- which((alleles$Pos == i) & (alleles$ID == recep))
           
      ARD <- alleles$Allele_Ratio[d]
      ARR <- alleles$Allele_Ratio[r]
      
      if (ARD > 0 & ARD < 1){
        Variance <- ARD*(1-ARD)/NT
        Log_estimate_sum <- Log_estimate_sum + dnorm(ARR, mean = ARD, sd = Variance, log = TRUE)
        print(Log_estimate_sum) 
      }
      
      
  }
}

BottleneckLikelihood("N4P7","N4P8",1)

#### Repeat function using the adjusted AR that excludes 0 and 1
BottleneckLikelihood_ad <- function(c,a,b){
  donor <- a
  recep <- b
  Log_estimate_sum <- 0
  NT <- c
  for (i in sites){
    d <- which((alleles_ad$Pos == i) & (alleles_ad$ID == donor))
    r <- which((alleles_ad$Pos == i) & (alleles_ad$ID == recep))
    
    ARD <- alleles_ad$Allele_Ratio_adj[d]
    ARR <- alleles_ad$Allele_Ratio_adj[r]
    

      Variance <- ARD*(1-ARD)/NT
      Log_estimate_sum <- Log_estimate_sum + dnorm(ARR, mean = ARD, sd = Variance, log = TRUE)
  }
  return(-Log_estimate_sum)
}

BottleneckLikelihood_ad(2,"N4P8","N4P9")

###################################################################################
# Infer optimal bottleneck size accounting for stochastic change in host (leonard et al)
### Convert the table to long format and plot the graph of likelihood vs NT

Stochastic_bottleneck <- function(c,a,b){
  donor <- a
  recep <- b
  Log_estimate_sum <- 0
  NT <- as.integer(c)
  range <- NT-1
  for (i in sites){
    d <- which((alleles_ad$Pos == i) & (alleles_ad$ID == donor))
    r <- which((alleles_ad$Pos == i) & (alleles_ad$ID == recep))
    
    ARD <- alleles_ad$Allele_Ratio_adj[d]
    ARR <- alleles_ad$Allele_Ratio_adj[r]
    if(0.1 < ARD & ARD < 0.9){
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

Stochastic_bottleneck_allsites <- function(c,a,b){
  donor <- a
  recep <- b
  Log_estimate_sum <- 0
  NT <- as.integer(c)
  range <- NT-1
  for (i in sites){
    d <- which((alleles_ad$Pos == i) & (alleles_ad$ID == donor))
    r <- which((alleles_ad$Pos == i) & (alleles_ad$ID == recep))
    
    ARD <- alleles_ad$Allele_Ratio_adj[d]
    ARR <- alleles_ad$Allele_Ratio_adj[r]
      for (k in 1:range){
        binom <- pbinom(k, size=NT, prob=ARD)
        betap <- pbeta(ARR, k, NT-k, ncp = 0, lower.tail = TRUE, log.p = FALSE)
        product <- binom * betap
        
        Log_estimate_sum <- Log_estimate_sum + product
        
        #print(c(i,NT,k,binom,betap,-Log_estimate_sum))
    }
  }
  return(-Log_estimate_sum)
}


BottleneckLikelihood_ad(200,"N4P8","N4P9")
Stochastic_bottleneck(200,"N4P8","N4P9")
optim(par=2,
      fn = Stochastic_bottleneck,
      a = "N4P8",
      b = "N4P9", lower = 2, upper = 10000, method = "Brent")


###################################################################################
### Create a loop that find the optimum NT for each transmission pair
# Step 1: subset the sNP distance for only transmisison steps

dists_trans <- dists %>% filter(timestep == 1)

for (i in 1:nrow(dists_trans)){
  dlevel <- dists_trans$ID1[i] 
  dID <- levels(dists_trans$ID1[i] )[dlevel]
  
  rlevel <- dists_trans$ID2[i] 
  rID <- levels(dists_trans$ID2[i] )[rlevel]
  
  result <- optim(par=1,
    fn = BottleneckLikelihood_ad,
    a = dID,
    b = rID, lower = 1, upper = 10000, method = "Brent")
  dists_trans$NTmax[i] <- result$par
  dists_trans$NTloglik[i] <- -result$value

Sresult <- optim(par=2,
                fn = Stochastic_bottleneck_allsites,
                a = dID,
                b = rID, lower = 1, upper = 10000, method = "Brent")
dists_trans$S_all_NTmax[i] <- Sresult$par
dists_trans$S_all_NTloglik[i] <- -Sresult$value
}

###################################################################################
### Create a loop that adds data on how many sites are mismatched and overall AR_change
# Step 1: Create the functions

Get_ARdiff<- function(a,b){
  #Initiate variables: count of matches, qruery indeces (qri)
  AR_diff <- numeric(0)
  qri1 <- ""
  qri2 <- ""
  qrid1 <- levels(a)[a]
  qrid2 <- levels(b)[b]
  for (qrpos in sites) {
    qri1 <- which(alleles_ad$Pos == qrpos & alleles_ad$ID == qrid1 )
    qri2 <- which(alleles_ad$Pos == qrpos & alleles_ad$ID == qrid2 )
    tmp <- alleles_ad$Allele_Ratio[qri1] - alleles_ad$Allele_Ratio[qri2]
    
    AR_diff <- c(AR_diff,abs(tmp))
  }
  sum(AR_diff)
}

Count_AR_changes <- function(a,b){
  #Initiate variables: count of matches, qmuery indeces (qmi)
  nmismatches <- 0
  qmi1 <- ""
  qmi2 <- ""
  qmid1 <- a
  qmid2 <- b
  for (qmpos in sites) {
    qmi1 <- which(alleles_ad$Pos == qmpos & alleles_ad$ID == qmid1 )
    qmi2 <- which(alleles_ad$Pos == qmpos & alleles_ad$ID == qmid2 )
    
    if (alleles_ad$Allele_Ratio[qmi1] != alleles_ad$Allele_Ratio[qmi2]) {
      nmismatches <- nmismatches + 1
    }
  }
  nmismatches
}


Get_mismatches <- function(a,b){
  #Initiate variables: count of matches, qmuery indeces (qmi)
  nmismatches <- 0
  qmi1 <- ""
  qmi2 <- ""
  qmid1 <- a
  qmid2 <- b
  for (qmpos in sites) {
    qmi1 <- which(alleles_ad$Pos == qmpos & alleles_ad$ID == qmid1 )
    qmi2 <- which(alleles_ad$Pos == qmpos & alleles_ad$ID == qmid2 )
    
    if (alleles_ad$SNP_type[qmi1] != alleles_ad$SNP_type[qmi2]) {
      nmismatches <- nmismatches + 1
    }
  }
  nmismatches
}


## Run a loop to get AR_Diff and number of mismatched states for all sites
for (i in 1:nrow(dists_trans)){
  dlevel <- dists_trans$ID1[i] 
  dID <- levels(dists_trans$ID1[i] )[dlevel]
  
  rlevel <- dists_trans$ID2[i] 
  rID <- levels(dists_trans$ID2[i] )[rlevel]
  
  #dists_trans$AR_diff[i] <- Get_ARdiff(dID,rID)
  #dists_trans$AR_changes[i] <- Count_AR_changes(dID,rID)
  dists_trans$State_changes[i] <- Get_mismatches(dID,rID)
}

########### Plot a box plot to see Nt distribution by AR sites changing

NT_Boxplot <- dists_trans %>% filter(AR_changes > 0) %>% 
  ggplot(aes(x = as.factor(AR_changes), y = NTmax)) + 
  geom_boxplot(aes(fill = AR_changes)) + geom_point() + 
  ylab("Optimum Bottleneck size") + xlab("Number of loci different") +
  theme_minimal() + theme(legend.position = "none")

dists_trans %>% filter(NTmax < 1000) %>% nrow()

NT_Histogram <- dists_trans %>% filter(AR_changes > 0) %>% 
  ggplot(aes(x = as.integer(NTmax))) + 
  xlab("Optimum Bottleneck size") + ylab("Probability Density") +
  geom_histogram(alpha = 0.2,fill=("blue"), col=("red"))

### Plot a panel figure with two using cowplot
plot_grid(NT_Boxplot, NT_Histogram)


### Plot graphs for NT from stochastic estimate
SNT_Boxplot <- dists_trans %>% filter(SNTmax < 1000) %>% 
  ggplot(aes(x = as.integer(SNTmax))) + 
  xlab("Optimum Bottleneck size") + ylab("Probability Density") +
  geom_histogram(alpha = 0.2,fill=("blue"), col=("red"))

SNT_Histogram <- dists_trans %>% filter(SNTmax < 1000) %>% 
  ggplot(aes(x = as.factor(AR_changes), y = as.integer(SNTmax))) + 
  geom_boxplot(aes(fill = AR_changes)) + geom_point() + 
  ylab("Optimum Bottleneck size") + xlab("Number of loci different") +
  theme_minimal() + theme(legend.position = "none")

plot_grid(SNT_Boxplot, SNT_Histogram)

dists_trans %>% filter(SNTmax < 1000) %>% nrow()

dists_trans %>% 
  ggplot(aes(x = as.integer(SNTmax))) + 
  xlab("Optimum Bottleneck size") + ylab("Probability Density") +
  geom_histogram(alpha = 0.2,fill=("blue"), col=("red"))


################### Old code and testing #############################3
### Convert the table to long format and plot the graph of likelihood vs NT

dists_long <- gather(dists_trans, Bottleneck, Loglikelihood, NT1:NT100, factor_key=TRUE)

dists_long <- dists_long %>% mutate(Bottleneck = str_replace(Bottleneck, "NT", ""))
#dists_long <- rename(dists_long, c("Bottleneck"="condition", "Loglikelihood"="measurement"))

max(dists_long %>% filter(ID1=="N1P4") %>% dplyr::select(Loglikelihood))

dists_long %>% filter(Ch1 == "N2") %>% ggplot(aes(x = as.numeric(Bottleneck), y = Loglikelihood)) + 
  geom_point() + facet_wrap(.~ID1, scales = "free") + xlab("Bottleneck size")

dists_long  %>% ggplot(aes(x = as.numeric(Bottleneck), y = Loglikelihood)) + 
  geom_boxplot() 

######## Create a loop to find max Bottleneck likelihood in each transmission
dists_long %>% filter(ID1 == "N1P13") %>%
  ggplot(aes(x = as.numeric(Bottleneck), y = Loglikelihood)) + 
  geom_point()

dists_trans %>% 
  ggplot(aes(x = as.integer(NTmax))) + 
  xlab("Optimum Bottleneck size") + ylab("Probability Density") +
  geom_histogram(alpha = 0.2,fill=("blue"), col=("red"))

