library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]  
percent <- ASV / mean(colSums(ASV)) *100

# nmds
percent.t <- t (percent)
nmds <-metaMDS(percent.t, trace=F, distance="bray", perm=100000)

# Make dataset
data <- cbind (nmds$points, DESIGN)

# ggplot
ggplot()+
geom_point(data=data, aes(x=MDS1,y=MDS2, color = Type, shape=Site,size=Season))+
scale_colour_manual(values = c("gray","black"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))

# Save
ggsave(file = "NMDS.png",width =5, height =4)