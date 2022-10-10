library(ggplot2)
library(stringr)
library(stringi)

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)
setwd("~/R/Analysis/4_Paddy/HNPF/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt",header=T)

# % table
ASV <- ASV.table [,1:(ncol(ASV.table)-7)] 
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)] 
percent <- ASV / mean(colSums(ASV)) *100
# Remove "f__" before phylum name
taxonomy <- data.frame(lapply(taxonomy, function(x){gsub(pattern="f__", replacement = "", x)}),stringsAsFactors = FALSE) 

# Aggregate
agrregated <- aggregate(percent, by=list(taxonomy$Family),FUN = sum,na.rm=F) 
row.names(agrregated)<-agrregated[,1]
agrregated <- agrregated[,-1]
agrregated <- data.frame(agrregated)
rowMeans <- rowMeans(agrregated)
agrregated <- cbind(agrregated,rowMeans)

# Main + <1% abund 
majors <- agrregated[agrregated[,"rowMeans"] > 1,] 
majors <- majors[order(majors$rowMeans,decreasing = T),]
minors <- agrregated[agrregated[,"rowMeans"] < 1,] 
Others <- colSums (minors)
selected <- majors
selected <- selected[,-ncol(selected)] 

# Make dataset
selected.t <- t (selected)
write.csv(selected.t, "aggregated.family.table.csv")
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

rownames<-rownames(selected)[nrow(selected):1]
data$Example <- factor(data$Example, levels = rownames)

# ggplot
ggplot (data,  mapping=aes(x=Type, y=Abundance, fill=Example))+
geom_bar(aes(), stat="identity", position="stack",color="black")+
theme_classic()+
theme(text=element_text(size=12,color="black"),axis.text=element_text(size=12,color="black"))+
labs (x="",y="Relative abundance (%)")+
facet_wrap(~Season) 

# Save
ggsave(file = "Rel.Abund.Family.png", width = 8,height=5)
