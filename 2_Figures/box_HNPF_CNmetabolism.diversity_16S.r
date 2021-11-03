library(vegan)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/16S")
setwd("~/R/Analysis/4_Paddy/HNPF/16S/Tax4Fun2/CNmetabolism")
C.table <- read.csv(file="C_metabolism.KOs.HNPF.csv",header=T,row.names=2)
N.table <- read.csv(file="N_metabolism.KOs.HNPF.csv",header=T,row.names=2)
C.table <- C.table[,-1]
N.table <- N.table[,-1]
C.table <- C.table[,-ncol(C.table)]
N.table <- N.table[,-ncol(N.table)]
C.table <- C.table[,-ncol(C.table)]
N.table <- N.table[,-ncol(N.table)]
setwd("~/R/Analysis/4_Paddy/HNPF/SSPN")

# For C metabolism
C.t <- t(C.table)
shannon <- diversity(C.t, index="shannon",base=2)
simpson <- diversity(C.t, index="simpson")
invsimpson <- diversity(C.t, index="invsimpson")
# fisher <- fisher.alpha(C.t)
# data<-cbind(shannon, simpson, invsimpson,fisher)
data<-cbind(shannon, simpson, invsimpson)
write.csv(data,"diversity.C.metabolism.csv")

data <- cbind(DESIGN,data)

ggplot(data) +
geom_boxplot(aes(y = shannon, x = Type, fill = Type))+   #Change
# scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="Shannon's diversity index",x="")+ # if you want to change the axis titles
facet_wrap(~Site)

ggsave(file = "Shannon.C.metabolism.png",width = 5, height = 4)


# For N metabolism
N.t <- t(N.table)
shannon <- diversity(N.t, index="shannon",base=2)
simpson <- diversity(N.t, index="simpson")
invsimpson <- diversity(N.t, index="invsimpson")
# fisher <- fisher.alpha(N.t)
# data<-cbind(shannon, simpson, invsimpson,fisher)
data<-cbind(shannon, simpson, invsimpson)
write.csv(data,"diversity.N.metabolism.csv")

data <- cbind(DESIGN,data)

ggplot(data) +
geom_boxplot(aes(y = shannon, x = Type, fill = Type))+   #Change
scale_fill_manual(values = c("white","black"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
scale_y_continuous(name = "Prokaryotic N-cycling \n gene diversity")+
facet_wrap(~Site)

ggsave(file = "Shannon.N.metabolism.png",width = 4, height = 2.5)
