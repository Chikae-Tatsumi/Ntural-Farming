library(vegan)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/16S/Tax4Fun2/CNmetabolism")
C.table <- read.csv(file="C_metabolism.KOs.HNPF.csv",header=T,row.names=2)
N.table <- read.csv(file="N_metabolism.KOs.HNPF.csv",header=T,row.names=2)
C.table <- C.table[,-1]
N.table <- N.table[,-1]

C.KOs <- C.table [,1:(ncol(C.table)-2)] 
colSums.C <- colSums(C.KOs)*100

N.KOs <- N.table [,1:(ncol(N.table)-2)] 
colSums.N <- colSums(N.KOs)*100

data <- cbind(DESIGN,colSums.C, colSums.N)

# C metabolism abudnance
ggplot(data) +
geom_boxplot(aes(y = colSums.C, x = Type, fill = Type))+   
scale_fill_manual(values = c("white","black"))+  
theme_classic()+
theme(text=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"))+
scale_y_continuous(name = "Prokaryotic C-cycling \n gene abundance (%)")+
facet_wrap(~Season)

ggsave(file = "total.C.metabolism.abund.png",width = 4, height = 3)

# N metabolism abudnance
ggplot(data) +
geom_boxplot(aes(y = colSums.N, x = Type, fill = Type))+   
scale_fill_manual(values = c("white","black"))+ 
theme_classic()+
theme(text=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"))+
scale_y_continuous(name = "Prokaryotic N-cycling \n gene abundance (%)")+
facet_wrap(~Season)

ggsave(file = "total.N.metabolism.abund.png",width = 4, height = 3)

# data save
data.save <- data[,-1:-ncol(DESIGN)]
colnames(data.save) <- c("C.metabolism.gege.rel.abund","N.metabolism.gege.rel.abund")
write.csv(data.save, "total.CNmetabolism.gene.abund.csv")
