#  see https://github.com/zdk123/SpiecEasi

# install SpiecEasi
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)  
library(Hmisc)  
library(Matrix)  

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF")
DESIGN <- read.csv("experimental_design.csv",header=T)

ASV.table <- read.table(file="16S/rarefied_ASV_table.txt",header=T)
ASV <- ASV.table [,1:(ncol(ASV.table)-6)]
taxonomy <- ASV.table [,(ncol(ASV.table)-5):ncol(ASV.table)]
percent <- ASV / mean(colSums(ASV)) *100
percent.t <- t(percent)

# To make minor phylum "Others"
taxonomy.minusp <- data.frame(lapply(taxonomy, function(x){gsub(pattern="p__", replacement = "", x)}),stringsAsFactors = FALSE)
rownames(taxonomy.minusp) <- rownames(taxonomy)
taxonomy <- taxonomy.minusp 
taxonomy$Phylum <- ifelse(is.na(taxonomy$Phylum), "unknown", taxonomy$Phylum) # replace NA with "unknown"
phylum <- aggregate(percent, by=list(taxonomy$Phylum),FUN = sum,na.rm=F) 
row.names(phylum)<-phylum[,1]
phylum <- phylum[,-1]
phylum <- data.frame(phylum)
rowMeans <- rowMeans(phylum) 
phylum <- cbind(phylum,rowMeans)
minor.phylum <- phylum[phylum[,"rowMeans"] < 2,] 
minor.phylum.list <- rownames(minor.phylum)
minor.phylum.list <- c(minor.phylum.list, "unknown")

for (i in 1:length (minor.phylum.list)){
taxonomy$Phylum <- gsub(minor.phylum.list[i],"Others",taxonomy$Phylum)}

sel.tax <- taxonomy

# Align taxonomy names
sel.tax$Phylum <- factor(sel.tax$Phylum)
others.n <- 1
others.n <- which(levels(sel.tax$Phylum)=="Others")
levels <- NULL
if (length(others.n) == 0){
     for (k in 1:(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}} else if (others.n == 1){
    for (k in 1:(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
    levels <- c(levels, "Others")} else{
for (k in 1:(others.n-1)){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
for (k in (others.n+1):(length(levels(sel.tax$Phylum)))){
    levels <- c(levels, levels(sel.tax$Phylum)[k])}
levels <- c(levels, "Others")
levels(sel.tax$Phylum) <- levels}

levels(sel.tax$Phylum)

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
sort <- sort(colMeans(subset),decreasing=T)
rank200 <- sort[200]
subset.t.filter <- subset[ ,colMeans(subset) >= rank200]
subset.t.filter <- subset.t.filter[,colSums(subset.t.filter)>0]
print(c(ncol(subset),"versus",ncol(subset.t.filter)))

# Calculate network (spiec.easi)
se <- spiec.easi(subset.t.filter, method='glasso', lambda.min.ratio=1e-2, nlambda=100, 
                          pulsar.params=list(rep.num=10)) 
saveRDS(se, paste("Network/SpiecEasi/top200/bacteria/",NAME[i],".se.rds", sep=""))

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
write.csv(data, paste("Network/SpiecEasi/top200/bacteria/",NAME[i],".metrics.csv", sep=""))
Results <- rbind (Results,data)
}
write.csv(Results, "Network/SpiecEasi/top200/bacteria/all.metrics.csv")
