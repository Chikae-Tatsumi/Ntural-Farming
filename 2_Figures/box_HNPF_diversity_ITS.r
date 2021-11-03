library(vegan)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/SSPN")

ASV <- ASV.table [,1:(ncol(ASV.table)-7)] # I have changed 7 into 6 because the table did not contain the Species column.
ASV.t <- t(ASV)
shannon <- diversity(ASV.t, index="shannon",base=2)
simpson <- diversity(ASV.t, index="simpson")
invsimpson <- diversity(ASV.t, index="invsimpson")
fisher <- fisher.alpha(ASV.t)
data<-cbind(shannon, simpson, invsimpson,fisher)
write.csv(data,"diversity.csv")

data <- cbind(DESIGN,shannon, simpson, invsimpson,fisher)

ggplot(data) +
geom_boxplot(aes(y = shannon, x = Type, fill = Type))+   #Change
scale_fill_manual(values = c("white","black"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
scale_y_continuous(name = "Fungal \n species diversity")+
facet_wrap(~Site)

ggsave(file = "Shannon.its.png",width = 4, height = 2.5)



ggplot(data) +
geom_boxplot(aes(y = fisher, x = Type, fill = Type))+   #Change
# scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Fisher's diversity index",x="")+ # if you want to change the axis titles
facet_wrap(~Site)

ggsave(file = "Fisher.png",width = 5, height = 4)