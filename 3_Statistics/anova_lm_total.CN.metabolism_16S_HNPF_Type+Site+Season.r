library (lme4)
library (lmerTest)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv")
setwd("~/R/Analysis/4_Paddy/HNPF/16S/Tax4Fun2")
METADATA <- read.csv("total.CNmetabolism.gene.abund.csv",header=T,row.names=1) #Change

design <- DESIGN
target <- METADATA

NAME1 <- "Type"
NAME2 <- "Site"
NAME3 <- "Season"

# ANOVA
Results <- NULL
for (i in 1:ncol(target)){
data <- cbind(design, as.numeric(target[,i]))
colnames (data)[ncol(data)] <- colnames(target)[i]
data <- na.omit (data)
anova <- anova(lm(scale(data[,ncol(data)])~ Type+Site+Season,data=data))
summary <- summary(lm(scale(data[,ncol(data)])~ Type+Site+Season,data=data))

estimate <- summary$coefficients[,1]

Fval <- anova[,4]
Pval <-anova[,5]
Bind <- c(estimate[2], Fval[1], Pval[1], estimate[3], Fval[2], Pval[2], estimate[4],Fval[3], Pval[3])
names(Bind) <- c(paste(NAME1,".estimate",sep=""),paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".estimate",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""),paste(NAME3,".estimate",sep=""),paste(NAME3,".F",sep=""),paste(NAME3,".P",sep=""))
Results <- rbind(Results, Bind)}

rownames(Results) <- colnames(target)

Results.asterisk <- Results
for (k in 1:ncol(Results)){
if(k%%3 == 0){
for (i in 1:nrow(Results)){
if (Results[i,k] < 0.001) {Results.asterisk[i,k] <- "***"
} else if (Results[i,k] < 0.01) {Results.asterisk[i,k] <- "**"
} else if (Results[i,k] < 0.05) {Results.asterisk[i,k] <- "*"
} else {Results.asterisk[i,k]<- "n.s"}}
}}

write.csv(Results.asterisk, "lm.total.CNmetabolism.withEstimate.Type+Site+Season.HNPF.csv") #Change

