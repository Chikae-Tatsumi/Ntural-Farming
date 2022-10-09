library(tidyverse)
library (dplyr)
library(ggplot2)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/ITS/FungalTrait")

fungaltrait <- read.csv(file="aggregated.fungaltrait.table.csv",header=T,row.names = 1)
data <- cbind (fungaltrait, DESIGN)

# Saprotroph abundance
ggplot(data) +
geom_boxplot(aes(y = Saprotroph, x = Type, fill = Type))+   
scale_fill_manual(values = c("white","black"))+  
theme_classic()+
theme(text=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"))+
scale_y_continuous(name = "Saprotroph fungal \n abundance (%)")+
facet_wrap(~Season)

ggsave("Saprotroph.png",width = 4, height = 3)
