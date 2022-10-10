library(vegan)
library(ggplot2)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

ASV <- ASV.table [,1:(ncol(ASV.table)-6)] 
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
simpson <- diversity(ASV.t, index="simpson")
invsimpson <- diversity(ASV.t, index="invsimpson")
data<-cbind(shannon, simpson, invsimpson)
write.csv(data,"diversity.csv")

data <- cbind(DESIGN,shannon, simpson, invsimpson)

ggplot(data) +
geom_boxplot(aes(y = shannon, x = Type, fill = Type))+   
scale_fill_manual(values = c("white","black"))+  
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
scale_y_continuous(name = "Prokaryotic \n species diversity")+
facet_wrap(~Season)

ggsave(file = "Shannon.16s.png",width = 4, height = 3)
