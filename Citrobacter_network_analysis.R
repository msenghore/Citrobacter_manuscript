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
library("ggraph")
library(igraph)
setwd("~/Projects/C._rodentium/gvcf_reanalysis/Final_analysis/")

dists <- read.csv("../Transmission_likelihood/Pairwise_comparisons.csv")
alleles <- read.csv("All_snps_summary_long_clean.csv")

Metadata <- alleles %>% dplyr::select("ID","Chain","Mouse_No") %>% unique()

### A function that takes a isolates, 
##### 1) ranks comparisons from isolates within 3 prior steps
##### 2) Takes the highest likelihood
##### 3) Breaks ties by time
##### 4) Suggests a likeliest transmission pair

for (i in 1:nrow(Metadata)){
  qID <- Metadata$ID[i]
  qM <- Metadata$Mouse_No[i]
  qMl <- qM - 3
  
  ###Only consider transmissions that are within potential donors within 4 steps
  ###### i.e. for mouse number 8, only isolates from mouses 5-8 will be considered
  tmp <- dists %>% filter(ID1 == qID | ID2 == qID, Mouse1 <= qM , Mouse2 <= qM, Mouse1 >= qMl, Mouse2 >= qMl)
  
  ### Rank all comparisons based on transmission likelihood
  tmp$rank <- rank(-tmp$P_trans_ARmean, ties.method = "min")
  
  ## Identify number of top hits
  ntops <- length(which(tmp$rank == 1))
  
  if(ntops == 1){
    qmatch <- which(tmp$rank == 1)
    
    if(tmp$ID1[qmatch] == qID){
      Metadata$Donor[i] <- tmp$ID2[qmatch]
    }else{Metadata$Donor[i] <- tmp$ID1[qmatch]}
  }
  else{Metadata$Donor[i] <- ntops}
  
}

Network_table <- Metadata %>% dplyr::select("Donor","ID") %>% filter(Donor != "multiple")

Net_meta <- Metadata %>% filter(ID %in% c(Network_table$Donor,Network_table$ID))
##3 plot network

Net_basic <- graph_from_data_frame(d=Network_table, directed=T) 
plot(Net_basic, vertex.size = 5)

# Make a palette of 3 colors
library(RColorBrewer)
coul  <- brewer.pal(4, "Set1") 

# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(networkST8)$Platoon))]
V(networkST8)$shape <- ifelse(V(networkST8)$Company == "A150", "square", "circle")
as.factor(V(networkST8)$Company)

c2 = cluster_leading_eigen(networkST8) 
coords = layout_with_fr(networkST8)
# Make the plot
#plot(c2, networkST8, layout=coords , vertex.color=my_color, vertex.size = 7, vertex.shapes = my_shape, vertex.label = V(networkST8)$Subject, edge.weight = E(networkST8)$Dist)
plot(networkST8, vertex.color=my_color, vertex.size = 7, vertex.shapes = V(networkST8)$shape, vertex.label = V(networkST8)$patient, edge.weight = E(networkST8)$Dist, vertex.label.dist = 1)
legend("bottomleft", legend=levels(as.factor(V(networkST8)$Platoon)) , title = "Platoon" , col = coul , bty = "n", pch=20 , pt.cex = 2.5, cex = 1, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
legend("bottomleft", legend=c("A150","C219"), title = "Company" ,bty = "n", pch=c(15,16) , pt.cex = 2.5, cex = 1, text.col="black" , horiz = FALSE, inset = c(0, 0.1))


components_25 <- components(network25, mode = c("weak", "strong"))



