#  see https://github.com/zdk123/SpiecEasi

# install SpiecEasi
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)  

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)

# percent
ASV.table <- read.table(file="ITS/rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-7)]
taxonomy <- ASV.table [,(ncol(ASV.table)-6):ncol(ASV.table)]
percent <- ASV / mean(colSums(ASV)) *100
percent.t.ITS <- t(percent)
colnames(percent.t.ITS) <- paste("Fungi_", colnames(percent.t.ITS), sep="")

# To make minor phylum "Others"
taxonomy.minusp <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "Fungi-", x)}),stringsAsFactors = FALSE)
rownames(taxonomy.minusp) <- rownames(taxonomy)
taxonomy <- taxonomy.minusp 
taxonomy$Phylum <- ifelse(is.na(taxonomy$Phylum), "unknown", taxonomy$Phylum) # replace NA wwith "unknown"
phylum <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F) 
row.names(phylum)<-phylum[,1]
phylum <- phylum[,-1]
phylum <- data.frame(phylum)
rowMeans <- rowMeans(phylum) 
phylum <- cbind(phylum,rowMeans)
minor.phylum <- phylum[phylum[,"rowMeans"] < 10,] 
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
sel.tax.ITS <- sel.tax.ITS[,-ncol(sel.tax.ITS)]
rownames(sel.tax.ITS) <- paste("Fungi_", rownames(sel.tax.ITS), sep="")

# percent
ASV.table <- read.table(file="16S/rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]
percent <- ASV / mean(colSums(ASV)) *100
percent.t.16S <- t(percent)
colnames(percent.t.16S) <- paste("Prokaryote_", colnames(percent.t.16S), sep="")

percent.t <- cbind (percent.t.16S, percent.t.ITS)

# To make minor phylum "Others"
taxonomy.minusp <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE)
rownames(taxonomy.minusp) <- rownames(taxonomy)
taxonomy <- taxonomy.minusp 
taxonomy$Phylum <- paste("Prokaryote-",taxonomy$Phylum,sep="")
taxonomy$Phylum <- ifelse(is.na(taxonomy$Phylum), "unknown", taxonomy$Phylum) # replace NA with "unknown"
phylum <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F) 
row.names(phylum)<-phylum[,1]
phylum <- phylum[,-1]
phylum <- data.frame(phylum)
rowMeans <- rowMeans(phylum) 
phylum <- cbind(phylum,rowMeans)
minor.phylum <- phylum[phylum[,"rowMeans"] < 10,] 
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
rownames(sel.tax.16S) <- paste("Prokaryote_", rownames(sel.tax.16S), sep="")

sel.tax <- rbind(sel.tax.16S, sel.tax.ITS)

# Subset
ag.perc <- aggregate(percent.t, by=list(DESIGN$Site,DESIGN$Type),FUN = sum,na.rm=F)
names <- ag.perc [,1:2]
colnames(names) <- c("Site","Type") 

bind <- cbind (percent.t,DESIGN)
subset.list <- list()
subset.list[[1]] <- subset(percent.t, bind$Type=="CF"&bind$Site=="Site_A") 
subset.list[[2]] <- subset(percent.t, bind$Type=="CF"&bind$Site=="Site_B") 
subset.list[[3]] <- subset(percent.t, bind$Type=="CF"&bind$Site=="Site_C") 
subset.list[[4]] <- subset(percent.t, bind$Type=="NF"&bind$Site=="Site_A") 
subset.list[[5]] <- subset(percent.t, bind$Type=="NF"&bind$Site=="Site_B") 
subset.list[[6]] <- subset(percent.t, bind$Type=="NF"&bind$Site=="Site_C") 

NAME <- c("SiteA.CF","SiteB.CF","SiteC.CF","SiteA.NF","SiteB.NF","SiteC.NF")

Results <- NULL
for (i in 1:6){
subset <- subset.list[[i]]

# Filter to pick up parts of ASVs
subset.ITS <- subset [,(ncol(percent.t.16S)+1):ncol(percent.t)]
subset.16S <- subset  [,1:ncol(percent.t.16S)]

sort.ITS <- sort(colMeans(subset.ITS),decreasing=T)
rank200.ITS <- sort.ITS[200]
subset.t.filter.ITS <- subset.ITS[ ,colMeans(subset.ITS) >= rank200.ITS] 
subset.t.filter.ITS <- subset.t.filter.ITS[,colSums(subset.t.filter.ITS)>0]

sort.16S <- sort(colMeans(subset.16S),decreasing=T)
rank200.16S <- sort.16S[200]
subset.t.filter.16S <- subset.16S[ ,colMeans(subset.16S) >= rank200.16S] 
subset.t.filter.16S <- subset.t.filter.16S[,colSums(subset.t.filter.16S)>0]

subset.t.filter <- cbind (subset.t.filter.ITS, subset.t.filter.16S)

subset.filter <- cbind(colnames(subset.t.filter), colMeans(subset.t.filter))
subset.filter <- data.frame(subset.filter)
colnames(subset.filter) <- c("ID","colMeans")
sel.tax$ID <- rownames(sel.tax)
tax.info <- merge(subset.filter, sel.tax, by="ID",sort=F, all.x=TRUE)
rownames(tax.info) <- tax.info$ID
saveRDS(tax.info, paste("Network/SpiecEasi/top200/fungi_bacteria/",NAME[i],".tax.info.rds", sep=""))
print(c(ncol(subset),"versus",ncol(subset.t.filter)))

# Calculate network (spiec.easi)
se <- spiec.easi(subset.t.filter, method='glasso', lambda.min.ratio=1e-2, nlambda=100, 
pulsar.params=list(rep.num=10))
saveRDS(se, paste("Network/SpiecEasi/top200/fungi_bacteria/",NAME[i],".se.rds", sep=""))

ig.gl     <- adj2igraph(se$refit$stars)
getStability(se) # This should be <0.05. https://github.com/zdk123/SpiecEasi

# Calculate Network metrics
deg <- degree(ig.gl,mode="all")
deg.max <- max(deg)
pat <- average.path.length(ig.gl, unconnected = F)
tra <- transitivity(ig.gl, type = "average")
cohesion <- vertex_connectivity(ig.gl)
gsize <- gsize(ig.gl)
edge_density <- round(edge_density(ig.gl),digit=5)
betweenness.max <- max(betweenness(ig.gl, directed=T, weights=NA))
deg.ave <- mean(deg)
deg.10 <- length(subset (deg, deg>10))
data <- cbind(edge_density, gsize, cohesion,tra,pat,deg.max, betweenness.max, deg.ave, deg.10)
write.csv(data, paste("Network/SpiecEasi/top200/fungi_bacteria/",NAME[i],".metrics.csv", sep=""))
Results <- rbind (Results,data)
}
write.csv(Results, "Network/SpiecEasi/top200/fungi_bacteria/all.metrics.csv")

