library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
METADATA <- read.csv(file = "metadata.csv",header=T)

setwd("~/R/Analysis/4_Paddy/HNPF/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]  
percent <- ASV / mean(colSums(ASV)) *100

# nmds
percent.t <- t (percent)
nmds <-metaMDS(percent.t, trace=F, distance="bray", perm=100000)
env <- cbind(METADATA [,1:4], METADATA[,6])
colnames(env) <- c("C", "N", "CN", "Moisture","pH")
env.fit <- envfit(nmds, env, perm=100000, na.rm=TRUE)

# Make dataset
data <- cbind (nmds$points, DESIGN)
env.values <- (env.fit$vectors[1]$arrows)
pval <- (env.fit$vectors[4]$pvals)
env.arrows <- cbind (env.values, pval)
env.arrows <- data.frame(subset (env.arrows, pval < 0.05))

# ggplot
rate = 1 
ggplot()+
geom_point(data=data, aes(x=MDS1,y=MDS2, color = Type, shape=Site,size=Season))+
scale_colour_manual(values = c("gray","black"))+  
geom_segment(data=env.arrows, aes(x = 0, y = 0, xend = (NMDS1*rate), yend = (NMDS2*rate)), arrow = arrow(length = unit(0.3,"cm")),color="black")+
geom_text_repel(data=env.arrows, aes(x=(NMDS1*rate), y=(NMDS2*rate), label=rownames(env.arrows)),  size=6, color="black") +
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))

# Save
ggsave(file = "NMDS_withArrows.png",width =5, height =4)
