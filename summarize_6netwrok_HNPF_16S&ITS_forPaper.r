library(igraph)  
library(Hmisc)  
library(Matrix)  

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)

# percent
ASV.table <- read.table(file="ITS/rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]
percent <- ASV / mean(colSums(ASV)) *100
percent.t.ITS <- t(percent)

# To make minor phylum "Others"
taxonomy.minusp <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "Fungi-", x)}),stringsAsFactors = FALSE)
rownames(taxonomy.minusp) <- rownames(taxonomy)
taxonomy <- taxonomy.minusp 
taxonomy$Phylum <- ifelse(is.na(taxonomy$Phylum), "unknown", taxonomy$Phylum) # replace NA into "unknown"
phylum <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F) 
row.names(phylum)<-phylum[,1]
phylum <- phylum[,-1]
phylum <- data.frame(phylum)
rowMeans <- rowMeans(phylum) 
phylum <- cbind(phylum,rowMeans)
minor.phylum <- phylum[phylum[,"rowMeans"] < 10,] # Change
minor.phylum.list <- rownames(minor.phylum)
minor.phylum.list <- c(minor.phylum.list, "unknown")

for (i in 1:length (minor.phylum.list)){
taxonomy$Phylum <- gsub(minor.phylum.list[i],"Fungi-Others",taxonomy$Phylum)}

sel.tax <- taxonomy

# Align taxonomy names
sel.tax$Phylum <- factor(sel.tax$Phylum)
others.n <- 1
others.n <- which(levels(sel.tax$Phylum)=="Fungi-Others")
levels <- NULL
if (length(others.n) == 0){
     for (k in 1:(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}} else if (others.n == 1){
    for (k in 1:(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
    levels <- c(levels, "Fungi-Others")} else{
for (k in 1:(others.n-1)){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
for (k in (others.n+1):(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
levels <- c(levels, "Fungi-Others")
levels(sel.tax$Phylum) <- levels}

levels(sel.tax$Phylum)
sel.tax.ITS <- sel.tax

# percent
ASV.table <- read.table(file="16S/rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]
percent <- ASV / mean(colSums(ASV)) *100
percent.t.16S <- t(percent)

percent.t <- cbind (percent.t.16S, percent.t.ITS)

# To make minor phylum "Others"
taxonomy.minusp <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE)
rownames(taxonomy.minusp) <- rownames(taxonomy)
taxonomy <- taxonomy.minusp 
taxonomy$Phylum <- paste("Prokaryote-",taxonomy$Phylum,sep="")
taxonomy$Phylum <- ifelse(is.na(taxonomy$Phylum), "unknown", taxonomy$Phylum) # replace NA into "unknown"
phylum <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F) 
row.names(phylum)<-phylum[,1]
phylum <- phylum[,-1]
phylum <- data.frame(phylum)
rowMeans <- rowMeans(phylum) 
phylum <- cbind(phylum,rowMeans)
minor.phylum <- phylum[phylum[,"rowMeans"] < 10,] # Change
minor.phylum.list <- rownames(minor.phylum)
minor.phylum.list <- c(minor.phylum.list, "unknown")

for (i in 1:length (minor.phylum.list)){
taxonomy$Phylum <- gsub(minor.phylum.list[i],"Prokaryote-Others",taxonomy$Phylum)}

sel.tax <- taxonomy

# Align taxonomy names
sel.tax$Phylum <- factor(sel.tax$Phylum)
others.n <- 1
others.n <- which(levels(sel.tax$Phylum)=="Prokaryote-Others")
levels <- NULL
if (length(others.n) == 0){
     for (k in 1:(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}} else if (others.n == 1){
    for (k in 1:(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
    levels <- c(levels, "Prokaryote-Others")} else{
for (k in 1:(others.n-1)){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
for (k in (others.n+1):(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
levels <- c(levels, "Prokaryote-Others")
levels(sel.tax$Phylum) <- levels}

levels(sel.tax$Phylum)
sel.tax.16S <- sel.tax

sel.tax <- rbind(sel.tax.16S, sel.tax.ITS[,-ncol(sel.tax.ITS)])


# Subset
ag.perc <- aggregate(percent.t, by=list(DESIGN$Site,DESIGN$Type),FUN = sum,na.rm=F)
names <- ag.perc [,1:2]
colnames(names) <- c("Site","Type") 

bind <- cbind (percent.t,DESIGN)
subset.list <- list()
subset.list[[1]] <- subset(percent.t, bind$Type=="CF"&bind$Site=="Site_A") # Change
subset.list[[2]] <- subset(percent.t, bind$Type=="CF"&bind$Site=="Site_B") # Change
subset.list[[3]] <- subset(percent.t, bind$Type=="CF"&bind$Site=="Site_C") # Change
subset.list[[4]] <- subset(percent.t, bind$Type=="NF"&bind$Site=="Site_A") # Change
subset.list[[5]] <- subset(percent.t, bind$Type=="NF"&bind$Site=="Site_B") # Change
subset.list[[6]] <- subset(percent.t, bind$Type=="NF"&bind$Site=="Site_C") # Change

NAME <- c("SiteA.CF","SiteB.CF","SiteC.CF","SiteA.NF","SiteB.NF","SiteC.NF")

Results <- NULL
for (i in 1:6){
subset <- subset.list[[i]]
# Filter to pick up parts of ASVs
subset.t.filter <- subset[ ,colMeans(subset) >= 0.05] # To pick up >0.05% ASVs
subset.t.filter <- subset.t.filter[,colSums(subset.t.filter)>0]
print(c(ncol(subset),"versus",ncol(subset.t.filter)))

# Calculate network
percent.cor <- rcorr(as.matrix(subset.t.filter), type="spearman")
percent.pval <- forceSymmetric(percent.cor$P) # Self-correlation as NA
# Select only the taxa for the filtered ASVs by using rownames of percent.pval
# sel.tax <- taxonomy[rownames(percent.pval),,drop=FALSE]
# Sanity check --> should be "[1] TRUE"
# all.equal(rownames(sel.tax), rownames(percent.pval))

p.yes <- percent.cor$P<0.05
r.yes <- percent.cor$r>0
r.high <- percent.cor$r>0.6
r.val <- percent.cor$r # select all the correlation values 
p.r.yes = p.yes*r.yes*r.val*r.high
adjm<-p.r.yes 

net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)

# hs <- hub_score(net.grph, weights=NA)$vector　
# as <- authority_score(net.grph, weights=NA)$vector
# pr <- page.rank(net.grph,directed=F)$vector　
deg <- degree(net.grph, mode="all")　
deg.max <- max(deg)
pat <- average.path.length(net.grph, unconnected = F)
# tra <- transitivity(net.grph, type = "global")
tra <- transitivity(net.grph, type = "average")
cohesion <- vertex_connectivity(net.grph)
gsize <- gsize(net.grph)
edge_density <- round(edge_density(net.grph),digit=5)
betweenness.max <- max(betweenness(net.grph, directed=T, weights=NA))
deg.ave <- mean(deg)
deg.10 <- length(subset (deg, deg>10))
data <- cbind(edge_density, gsize, cohesion,tra,pat,deg.max, betweenness.max, deg.ave, deg.10)
Results <- rbind (Results,data)

# Illustrate
col=rainbow(length(levels(sel.tax$Phylum)))
plot.igraph(net.grph, vertex.size=deg*0.1,vertex.label=NA, vertex.color=col[unclass(sel.tax$Phylum)],layout=layout.kamada.kawai)
# plot(net.grph, vertex.size=deg*0.15,vertex.label=NA,vertex.color=col[unclass(sel.tax$Phylum)],layout=layout.random)
# plot(net.grph, vertex.size=deg*0.15,vertex.label=NA,vertex.color=col[unclass(sel.tax$Phylum)],layout=layout.fruchterman.reingold)
if (i==1){
legend(x = -4, y = 2, legend = levels(sel.tax$Phylum), cex=3, pch = 19, col = col, bty = "n",  y.intersp = 0.05, x.intersp = 0.05)}
text(x=0,y=-0.4,paste("Edge: ", format(gsize,digit=3)),cex=3)
text(x=0,y=-0.5,paste("Density: ", format(edge_density,digit=3)),cex=3)
title(NAME[i],cex=2) #Change

# Save 
if (i==1){
dev.copy(png, file=paste("Network_16SITS/forPaper/Network.HNPF.16S.",NAME[i],".png",sep=""), height=800, width=1400)
dev.off()} else{
dev.copy(png, file=paste("Network_16SITS/forPaper/Network.HNPF.16S.",NAME[i],".png",sep=""), height=800, width=800)
dev.off()}}

Results <- cbind(names,Results)
colnames(Results) <- c("Site","Type","edge.density","edge.number", "cohesion","clustering.coefficient","ave.path.length","maximum.degree", "max.betweenness", "average.degree", ">10.degree")

write.csv(Results, "Network_16SITS/forPaper/network.indicates.season.composited.csv")
