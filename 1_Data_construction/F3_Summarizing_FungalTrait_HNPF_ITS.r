library(tidyverse)
library (dplyr)
library(ggplot2)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/ITS/FungalTrait")

fungaltrait.table <- read.csv(file="fungaltrait.table.csv",header=T,row.names = 1)
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100
guild$lifestyle <- paste(guild$primary_lifestyle, "_",guild$secondary_lifestyle)

# For ECM
aggregated <- aggregate(percent, by=list(guild$lifestyle),FUN = sum,na.rm=F) 
row.names(aggregated)<-aggregated[,1]
aggregated <- aggregated[,-1]
aggregated <- data.frame(aggregated)
rows <- grep ("ectomycorrhizal", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
ECM <- data.frame(rowSums (subset.t))
colnames (ECM) [1] <- "ECM"

# For Saprotroph
rows <- grep ("saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
Saprotroph <- data.frame(rowSums (subset.t))
colnames (Saprotroph) [1] <- "Saprotroph" 

# For Pathotroph
rows <- grep ("patho", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
Pathotroph <- data.frame(rowSums (subset.t))
colnames (Pathotroph) [1] <- "Pathotroph"

# Plant ptatogenic capacity
pat.cap <- cbind(percent, guild$Plant_pathogenic_capacity_template)
colnames(pat.cap)[ncol(pat.cap)] <- "Plant_pathogenic_capacity_template"
pat.cap <- na.omit(pat.cap)
pat.cap[pat.cap$Plant_pathogenic_capacity_template == "",] <-NA
pat.cap <- na.omit (pat.cap)
colSums <- colSums(pat.cap[,-ncol(pat.cap)])
Plant_pathogenic_capacity <- data.frame(colSums)
colnames (Plant_pathogenic_capacity) [1] <- "Plant_pathogenic_capacity"

# Animal_biotrophic_capacity
pat.cap <- cbind(percent, guild$Animal_biotrophic_capacity_template)
colnames(pat.cap)[ncol(pat.cap)] <- "Animal_biotrophic_capacity_template"
pat.cap <- na.omit(pat.cap)
pat.cap[pat.cap$Animal_biotrophic_capacity_template == "",] <-NA
pat.cap <- na.omit (pat.cap)
rows <- grep ("parasite", pat.cap$Animal_biotrophic_capacity_template)
subset <- pat.cap[rows,]
subset <- subset[,-ncol(subset)]
Animal_parasite <- data.frame(colSums (subset))
colnames (Animal_parasite) [1] <- "Animal_parasite"

# For soil saprotroph
rows <- grep ("soil_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
soil_saprotroph <- data.frame(rowSums (subset.t))
colnames (soil_saprotroph) [1] <- "soil_saprotroph"

# For wood saprotroph
rows <- grep ("wood_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
wood_saprotroph <- data.frame(rowSums (subset.t))
colnames (wood_saprotroph) [1] <- "wood_saprotroph"

# For litter saprotroph
rows <- grep ("litter_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
litter_saprotroph <- data.frame(rowSums (subset.t))
colnames (litter_saprotroph) [1] <- "litter_saprotroph"

# For dung saprotroph
rows <- grep ("dung_saprotroph", rownames (aggregated)) 
subset <- aggregated[rows,]
subset.t <- t(subset)
dung_saprotroph <- data.frame(rowSums (subset.t))
colnames (dung_saprotroph) [1] <- "dung_saprotroph"

# Summarize
Result <- cbind (ECM, Saprotroph, Pathotroph,Plant_pathogenic_capacity, Animal_parasite,soil_saprotroph, wood_saprotroph, litter_saprotroph, dung_saprotroph)
write.csv(Result, "aggregated.fungaltrait.table.csv")
