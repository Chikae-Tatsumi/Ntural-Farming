library(vegan)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/16S")

ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
ASV.t <- t(ASV)

# PerMANOVA
adonis <- adonis(ASV.t ~ Type*Season, data=DESIGN, permutations=10000
, strata = DESIGN$Site
) # CHANGE ME

# Save
Fval <- adonis[[1]][,4]
R2 <- adonis[[1]][,5]
Pval <-adonis[[1]][,6]
bind <- cbind(Fval, R2, Pval)
rownames (bind) <- rownames(adonis[[1]])
write.csv (bind, "PerMANOVA.csv")
