library(ggplot2)
library(vegan)
library(ggrepel)
library(qgraph)
library(psych)
library(GGally) 
library(corpcor)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
design <- read.csv(file = "experimental_design.csv",header=T)

setwd("~/R/Analysis/4_Paddy/HNPF/16S/Tax4Fun2")
KOs <- read.csv("functional_prediction.csv",header=T, row.names=1)
pathway <- read.csv("pathway_prediction.csv",header=T)

# aggregate KOs
KOs <- KOs[,-ncol(KOs)]
KOs.t <- t(KOs)
KOs.t <- data.frame(KOs.t)
pathway.values <- pathway[,2:(ncol(pathway)-3)]

ag.KOs <- aggregate(KOs.t, by=list(DESIGN$Type,DESIGN$Site, DESIGN$Season),FUN = mean,na.rm=F)
colnames(ag.KOs)[1:3] <- c("Type","Site","Season")
write.csv (ag.KOs,"aggregated.Tax4Fun2.KOs.csv")

# Count functional genes based on EC number
bg.bind <- cbind(KOs.t$K01188,KOs.t$K05349,KOs.t$K05350)
bx.bind <- cbind(KOs.t$K01198,KOs.t$K15920,KOs.t$K22268)
ag.bind <- cbind(KOs.t$K01187,KOs.t$K12047,KOs.t$K12316,KOs.t$K12317)
cel.bind <- cbind(KOs.t$K01179,KOs.t$K19357,KOs.t$K20542)
nag.bind <- cbind(KOs.t$K01207,KOs.t$K12373,KOs.t$K14459,KOs.t$K20730)
chi.bind <- cbind(KOs.t$K01183, KOs.t$K13381,KOs.t$K20547)
lap.bind <- cbind(KOs.t$K01255, KOs.t$K11142)
ure.bind <- cbind(KOs.t$K01427,KOs.t$K01428,KOs.t$K01429,KOs.t$K01430,KOs.t$K14048)
fix.bind <- cbind(KOs.t$K00531,KOs.t$K02586, KOs.t$K02591)
amo.bind <- KOs.t$K10944 # EC 1.14.99.39 amoA
nir.bind <- cbind(KOs.t$K00368, KOs.t$K15864)	# EC 1.7.2.1 nirK, nirS
nos.bind <- KOs.t$K00376 # EC 1.7.2.4 nosZ
ap.bind <- cbind(KOs.t$K01078, KOs.t$K01093,KOs.t$K03788,KOs.t$K09474,KOs.t$K14379,KOs.t$K14394,KOs.t$K14395,KOs.t$K14410,KOs.t$K19283,KOs.t$K19284)
cbh.bind <- cbind(KOs.t$K01225, KOs.t$K19668)
Arginase <- KOs.t$K01476
PhOx <- KOs.t$K05909
per.bind <- cbind(KOs.t$K00430, KOs.t$K10788,KOs.t$K11188,KOs.t$K12550,KOs.t$K19511,KOs.t$K21201)

BG <- rowSums(bg.bind)
BX <- rowSums(bx.bind)
AG <- rowSums(ag.bind)
Cellulase <- rowSums(cel.bind)
NAG <- rowSums(nag.bind)
Chitinase <- rowSums(chi.bind)
LAP <- rowSums(lap.bind)
Urease <- rowSums(ure.bind)
N.fix <-  rowSums(fix.bind)
amoA <- amo.bind
nirSK <-  rowSums(nir.bind)
nosZ <- nos.bind
aP <- rowSums(ap.bind)
CBH <- rowSums(cbh.bind)
Perox <- rowSums(per.bind)

data <- cbind(BG, AG, BX, Cellulase, CBH,Chitinase,NAG,LAP, Arginase,Urease, N.fix,amoA, nirSK, nosZ, aP)
if (length (PhOx)>0){data <- cbind(data,PhOx)}
if (length (Perox)>0){data <- cbind(data,Perox)}
rownames(data) <- colnames(KOs)
write.csv(data, "Tax4Fun2.CNgenes.csv")

# aggregate pathways
pathway.level3 <- aggregate(pathway.values, by=list(pathway$level3),FUN = sum,na.rm=F)
pathway.level2 <- aggregate(pathway.values, by=list(pathway$level3,pathway$level2),FUN = sum,na.rm=F)
pathway.level1 <- aggregate(pathway.values, by=list(pathway$level3,pathway$level2,pathway$level1),FUN = sum,na.rm=F)
colnames(pathway.level3)[1] <- "level3"
colnames(pathway.level2)[1:2] <- c("level3","level2")
colnames(pathway.level2)[1:3] <- c("level3","level2","level1")
pathway.level2 <- pathway.level2[order(pathway.level2$level3),]
write.csv (pathway.level3,"level3.Tax4Fun2.pathway.csv")
write.csv (pathway.level2 ,"level2.Tax4Fun2.pathway.csv")
write.csv (pathway.level1 ,"level1.Tax4Fun2.pathway.csv")

pathway.level2 <- pathway.level2[,-1]
pathway.level1 <- pathway.level1[,-1:-2]
pathway.level3.t <- t(pathway.level3[,-1])
pathway.level2.t <- t(pathway.level2[,-1])
pathway.level1.t <- t(pathway.level1[,-1])
colnames(pathway.level1.t) <- pathway.level1[,1]
colnames(pathway.level2.t) <- pathway.level2[,1]
colnames(pathway.level3.t) <- pathway.level3[,1]
ag.level3 <- aggregate(pathway.level3.t, by=list(DESIGN$Type,DESIGN$Site, DESIGN$Season),FUN = mean,na.rm=F)
ag.level2 <- aggregate(pathway.level2.t, by=list(DESIGN$Type,DESIGN$Site, DESIGN$Season),FUN = mean,na.rm=F)
ag.level1 <- aggregate(pathway.level1.t, by=list(DESIGN$Type,DESIGN$Site, DESIGN$Season),FUN = mean,na.rm=F)
colnames(ag.level1)[1:3] <- c("Type","Site","Season")
colnames(ag.level2)[1:3] <- c("Type","Site","Season")
colnames(ag.level3)[1:3] <- c("Type","Site","Season")

write.csv (ag.level3,"aggregated.level3.Tax4Fun2.pathway.csv")
write.csv (ag.level2 ,"aggregated.level2.Tax4Fun2.pathway.csv")
write.csv (ag.level1 ,"aggregated.level1.Tax4Fun2.pathway.csv")

