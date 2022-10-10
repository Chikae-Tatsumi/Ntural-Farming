library (lme4)
library (lmerTest)

setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv")
setwd("~/R/Analysis/4_Paddy/HNPF/16S/Tax4Fun2/CNgenes")
METADATA <- read.csv("Tax4Fun2.CNgenes.csv",header=T,row.names=1) 

design <- DESIGN
target <- METADATA

NAME1 <- "Type"
NAME2 <- "Season"
NAME3 <- "Site"

# ANOVA
Results <- NULL
for (i in 1:ncol(target)){
data <- cbind(design, as.numeric(target[,i]))
colnames (data)[ncol(data)] <- colnames(target)[i]
data <- na.omit (data)

Model <- lmer(scale(data[,ncol(data)])~ Type*Season+(1|Site),data=data)
anova <- anova(Model)
summary <- summary(Model)

estimate <- summary$coefficients[,1]

Fval <- anova[,5]
Pval <-anova[,6]
Bind <- c( estimate[2],Fval[1], Pval[1], Fval[2], Pval[2], Fval[3], Pval[3])
names(Bind) <- c("TypeNF.estimate", paste(NAME1,".F",sep=""),paste(NAME1,".P",sep=""),paste(NAME2,".F",sep=""),paste(NAME2,".P",sep=""),paste(NAME1,"*",NAME2,".F",sep=""),paste(NAME1,"*",NAME2,".P",sep=""))
Results <- rbind(Results, Bind)}

rownames(Results) <- colnames(target)

Results.asterisk <- Results
for (k in 3:ncol(Results)){
if(k%%2 == 1){
for (i in 1:nrow(Results)){
if (Results[i,k] < 0.001) {Results.asterisk[i,k] <- "***"
} else if (Results[i,k] < 0.01) {Results.asterisk[i,k] <- "**"
} else if (Results[i,k] < 0.05) {Results.asterisk[i,k] <- "*"
} else {Results.asterisk[i,k]<- "n.s"}}
}}

write.csv(Results.asterisk, "lmer.CNgenes.Type×Season+(Site).HNPF.csv") 

