library(vegan)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv(file = "experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/ITS")

ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
ASV.t <- t(ASV)

# If you want to remove NA data
# bind <- cbind (ASV.t, DESIGN)
# bind.rm <- na.omit (bind)
# ASV.t <- bind.rm [, 1:nrow(ASV.table)]
# DESIGN.rm <- na.omit (DESIGN)

# PerMANOVA
adonis <- adonis(ASV.t ~ Type+Site+Season, data=DESIGN, permutations=10000
# , strata = DESIGN$Group
) # CHANGE ME
# You can use DESIGN.rm instead of DESIGN if you removed NA data

# Save
Fval <- adonis[[1]][,4]
R2 <- adonis[[1]][,5]
Pval <-adonis[[1]][,6]
bind <- cbind(Fval, R2, Pval)
rownames (bind) <- rownames(adonis[[1]])
write.csv (bind, "PerMANOVA.csv")


# PerMANOVA
adonis <- adonis(ASV.t ~ Type*Site+Season, data=DESIGN, permutations=10000
# , strata = DESIGN$Group
) # CHANGE ME
# You can use DESIGN.rm instead of DESIGN if you removed NA data

# Save
Fval <- adonis[[1]][,4]
R2 <- adonis[[1]][,5]
Pval <-adonis[[1]][,6]
bind <- cbind(Fval, R2, Pval)
rownames (bind) <- rownames(adonis[[1]])
write.csv (bind, "PerMANOVA_multiple.csv")