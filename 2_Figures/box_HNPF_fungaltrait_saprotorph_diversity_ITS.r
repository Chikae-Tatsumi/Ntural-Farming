library(vegan)
library(ggplot2)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/ITS/FungalTrait")

fungaltrait.table <- read.csv(file="fungaltrait.table.csv",header=T,row.names = 1)
ASV <- fungaltrait.table [,1:(ncol(fungaltrait.table)-31)] 
guild <- fungaltrait.table [,(ncol(fungaltrait.table)-24):ncol(fungaltrait.table)] 
percent <- ASV / mean(colSums(ASV)) *100
guild$lifestyle <- paste(guild$primary_lifestyle, "_",guild$secondary_lifestyle)
setwd("~/R/Analysis/4_Paddy/HNPF/SSPN")

percent.table <- cbind (percent, guild)
sap.table<- percent.table[grep("saprotroph", guild$lifestyle),]

ASV.sap <- sap.table [,1:(ncol(sap.table)-26)] 
ASV.sap <- as.matrix(ASV.sap)
ASV.t <- t(ASV.sap)
shannon <- diversity(ASV.t, index="shannon",base=2)
simpson <- diversity(ASV.t, index="simpson")
invsimpson <- diversity(ASV.t, index="invsimpson")
# fisher <- fisher.alpha(ASV.t)
bind <- cbind(shannon, simpson, invsimpson)
write.csv(bind,"fungaltrait.sap.diversity.csv")

data <- cbind(DESIGN,bind)

# Saprotroph
ggplot(data) +
geom_boxplot(aes(y = shannon, x = Type, fill = Type))+   #Change
scale_fill_manual(values = c("white","black"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
scale_y_continuous(name = "Saprotroph fungal \n diversity")+
facet_wrap(~Site)

ggsave(file = "fungaltrait.saprotroph.diversity.png",width = 4, height = 3)


