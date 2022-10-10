library(ggplot2)
library(stringr)
library(stringi)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/16S")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]  
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]  
percent <- ASV / mean(colSums(ASV)) *100
# Remove "k__", "p__" and "c__" before phylum name
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="k__", replacement = "", x)}),stringsAsFactors = FALSE)
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE)
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="c__", replacement = "", x)}),stringsAsFactors = FALSE)

# Separate Bacteria and Archaea
percent.table <- cbind (percent, taxonomy)
bacteria.table <- subset (percent.table, percent.table$Kingdom == "Bacteria")
archaea.table <- subset (percent.table, percent.table$Kingdom == "Archaea")
bacteria.ASV <- bacteria.table [,1:(ncol(ASV.table)-6)] 
archaea.ASV <- archaea.table [,1:(ncol(ASV.table)-6)] 
bacteria.taxonomy <- bacteria.table  [,(ncol(ASV.table)-5):ncol(ASV.table)] 
archaea.taxonomy <- archaea.table  [,(ncol(ASV.table)-5):ncol(ASV.table)]  

# >1% bacterial phylum table
phylum.ag <- aggregate(bacteria.ASV, by=list(bacteria.taxonomy$Phylum),FUN = sum,na.rm=F)
row.names(phylum.ag)<- phylum.ag[,1]
phylum.ag <- phylum.ag[,-1]
phylum.ag <- data.frame(phylum.ag)
rowMeans <- rowMeans(phylum.ag)
phylum.ag <- cbind(phylum.ag,rowMeans)
major.phylum <- phylum.ag[phylum.ag[,"rowMeans"] > 1,]
major.phylum <- major.phylum[order(major.phylum$rowMeans,decreasing = T),]
others <- colSums(bacteria.ASV)-colSums (major.phylum[,-ncol(major.phylum)])
selected.phylum <- rbind (major.phylum, others) 
rownames (selected.phylum)[nrow(selected.phylum)] <- "Other_Bacteria"
selected.phylum <- selected.phylum[,-ncol(selected.phylum)] 
selected.phylum.t <- t (selected.phylum)
write.csv(selected.phylum.t, "aggregated.phylum.table.csv")

# 1%> archaeal table
a.phylum.ag <- aggregate(archaea.ASV, by=list(archaea.taxonomy$Phylum),FUN = sum,na.rm=F)
row.names(a.phylum.ag)<- a.phylum.ag[,1]
a.phylum.ag <- a.phylum.ag[,-1]
a.phylum.ag <- data.frame(a.phylum.ag)
rowMeans <- rowMeans(a.phylum.ag)
a.phylum.ag <- cbind(a.phylum.ag,rowMeans)
a.major.phylum <- a.phylum.ag[a.phylum.ag[,"rowMeans"] > 1,]
a.major.phylum <- a.major.phylum[order(a.major.phylum$rowMeans,decreasing = T),]
a.others <- colSums(archaea.ASV)-colSums (a.major.phylum[,-ncol(a.major.phylum)])
a.selected.phylum <- rbind (a.major.phylum, a.others) 
rownames (a.selected.phylum)[nrow(a.selected.phylum)] <- "Other_Archaea"
colnames (a.selected.phylum) <- colnames (a.phylum.ag[,-ncol(a.phylum.ag)])
a.selected.phylum <- a.selected.phylum[,-ncol(a.selected.phylum)]

# >1% class in >10% class tables
more.major.phylum <- phylum.ag[phylum.ag[,"rowMeans"] > 10,]
more.major.phylum <- more.major.phylum[order(more.major.phylum$rowMeans,decreasing = T),]
percent.table <- cbind (percent, taxonomy)
class.list <- list()
for (i in 1:nrow(more.major.phylum)){
subset <-subset(percent.table,percent.table$Phylum==rownames(more.major.phylum)[i])
class.list <- c(class.list, list(subset))}

selected.class.list <- list()
for (i in 1:length(class.list)){
class.ASV <- class.list[[i]][,1:(ncol(ASV.table)-6)] 
class.taxonomy <- class.list[[i]][,(ncol(ASV.table)-5):ncol(ASV.table)]  
class.ag <- aggregate(class.ASV, by=list(class.taxonomy$Class),FUN = sum,na.rm=F)
row.names(class.ag)<- class.ag[,1]
class.ag <- class.ag[,-1]
class.ag <- data.frame(class.ag)
rowMeans <- rowMeans(class.ag)
class.ag <- cbind(class.ag,rowMeans)
major.class <- class.ag[class.ag[,"rowMeans"] > 1,]
major.class <- major.class[order(major.class$rowMeans,decreasing = T),]
other.class <- colSums(class.ASV)-colSums (major.class[,-ncol(major.class)])
selected.class <- rbind (major.class, other.class) 
rownames (selected.class)[nrow(selected.class)] <- paste ("Other", rownames(more.major.phylum)[i],sep="_")
selected.class <- selected.class[,-ncol(selected.class)]
selected.class.list <- c(selected.class.list, list(selected.class))}

# Make dataset
phylum.class.table <- data.frame()
for (i in 1:nrow(more.major.phylum)){
phylum.class.table <- rbind(phylum.class.table,selected.class.list[[i]])}
phylum.class.table <- rbind (phylum.class.table, selected.phylum [((nrow(more.major.phylum)+1):nrow(selected.phylum)),])
phylum.class.table <- rbind (phylum.class.table, a.selected.phylum)
selected.t <- t (phylum.class.table)
write.csv(selected.t, "aggregated.phylum.class.table.csv")

bind <- cbind (selected.t, DESIGN)
category <- paste (bind$Type, bind$Season,sep="&") 
bind <- cbind (bind, category)
table <- aggregate(selected.t, by=list(bind$category),FUN = mean) 
split <- str_split (table[,1], "&",simplify = TRUE)
colnames(split) <- c("Type","Season") 
table <- cbind (split,table[2:(ncol(table))])
Type <- as.character(table[,1]) 
Season <- as.character(table[,2]) 
table <- data.frame(table)
data <- data.frame()
for (i in 3:(ncol(table))){
Abundance <-table[,i]
Example <- colnames (table)[i]
bind <- cbind(Type, Season, Abundance, Example) 
data<-rbind (data,bind)}
data$Abundance <- as.numeric(as.character(data$Abundance))

rownames<-colnames(table)[ncol(table):3]
data$Example <- factor(data$Example, levels = rownames)

# ggplot
ggplot (data,  mapping=aes(x=Type, y=Abundance, fill=Example))+
geom_bar(aes(), stat="identity", position="stack",color="black")+
scale_fill_manual(values = c("black","#454545","#919191"
,"lightgray","white","#dbffef","#ffe99e","#efdbff","#529eff","#c59900","#C77CFF","#85f1f9","#f9c785","#85f98d","#8d85f9","#faa49e"
,"#f0ffcb","#caff48","#a0e300","#7cb000","#466200","#223100"
,"#89fcff","#00bfc4", "#007478"
,"#fee9e7","#F8766D","#753c37"))+  
theme_classic()+
theme(text=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="",y="Abundance (%)")+
facet_wrap(~Season)

# Save
ggsave(file = "Rel.Abund.phylum.png", height=5,width = 8)