library(vegan)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/16S/Tax4Fun2/CNmetabolism")
C.table <- read.csv(file="C_metabolism.KOs.HNPF.csv",header=T,row.names=2)
N.table <- read.csv(file="N_metabolism.KOs.HNPF.csv",header=T,row.names=2)
C.table <- C.table[,-1]
N.table <- N.table[,-1]
setwd("~/R/Analysis/4_Paddy/HNPF/SSPN")

C.KOs <- C.table [,1:(ncol(C.table)-2)] # I have changed 7 into 6 because the table did not contain the Species column.
colSums.C <- colSums(C.KOs)*100

N.KOs <- N.table [,1:(ncol(N.table)-2)] # I have changed 7 into 6 because the table did not contain the Species column.
colSums.N <- colSums(N.KOs)*100

data <- cbind(DESIGN,colSums.C, colSums.N)

# C metabolism abudnance
ggplot(data) +
geom_boxplot(aes(y = colSums.C, x = Type, fill = Type))+   #Change
# scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=12,color="black"))+
labs (y="C metabolism gene abudnance (%)",x="")+ # if you want to change the axis titles
facet_wrap(~Site)

ggsave(file = "total.C.metabolism.abund.png",width = 5, height = 4)


# N metabolism abudnance
ggplot(data) +
geom_boxplot(aes(y = colSums.N, x = Type, fill = Type))+   #Change
scale_fill_manual(values = c("white","black"))+  # if you want to change the colors
theme_classic()+
theme(text=element_text(size=14,color="black"),axis.text=element_text(size=14,color="black"))+
scale_y_continuous(name = "Prokaryotic N-cycling \n gene abundance (%)")+
facet_wrap(~Site)

ggsave(file = "total.N.metabolism.abund.png",width = 4, height = 2.5)

# data save
data.save <- data[,-1:-ncol(DESIGN)]
colnames(data.save) <- c("C.metabolism.gege.rel.abund","N.metabolism.gege.rel.abund")
write.csv(data.save, "total.CNmetabolism.gene.abund.csv")
