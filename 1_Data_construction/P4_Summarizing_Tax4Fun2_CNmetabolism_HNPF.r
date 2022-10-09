library(ggplot2)
library(vegan)
library(ggrepel)

# Import files
C_metabolism <- read.csv(file = "~/R/Database/kegg_C_metabolism.csv",header=T)
N_metabolism <- read.csv(file = "~/R/Database/kegg_N_metabolism.csv",header=T)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/16S/Tax4Fun2")
METADATA <- read.csv("functional_prediction.csv",header=T) #Change

# merge
C.data <- merge( METADATA,C_metabolism, by ='KO')
N.data <- merge( METADATA,N_metabolism, by ='KO')

write.csv(C.data,"C_metabolism.KOs.HNPF.csv")
write.csv(N.data,"N_metabolism.KOs.HNPF.csv")