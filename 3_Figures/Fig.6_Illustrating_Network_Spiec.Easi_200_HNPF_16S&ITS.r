# install SpiecEasi
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)  

# Import files
setwd("~/R/Analysis/4_Paddy/HNPF/Network/SpiecEasi/top200/fungi_bacteria")
NAME <- c("SiteA.CF","SiteB.CF","SiteC.CF","SiteA.NF","SiteB.NF","SiteC.NF")

for (i in 1:6){
tax.info <- readRDS(paste(NAME[i],".tax.info.rds", sep=""))
se <- readRDS(paste(NAME[i],".se.rds", sep=""))
ig.gl     <- adj2igraph(se$refit$stars)
gsize <- gsize(ig.gl)
edge_density <- round(edge_density(ig.gl),digit=5)
deg <- degree(ig.gl,mode="all")

# Illustrate
col=rainbow(length(levels(tax.info$Phylum)))
plot.igraph(ig.gl, vertex.size=deg*0.3,vertex.label=NA, vertex.color=col[unclass(tax.info$Phylum)],layout=layout.kamada.kawai)
if (i==1){
legend(x = -5, y = 1, legend = levels(tax.info$Phylum), cex=3, pch = 19, col = col, bty = "n",  y.intersp = 0.1, x.intersp = 0.1)}
text(x=0,y=-0.5,paste("Edge: ", format(gsize,digit=3)),cex=4)
text(x=0,y=-0.7,paste("Density: ", format(edge_density,digit=3)),cex=4)
title(NAME[i],cex=2) #Change

# Save 
if (i==1){
dev.copy(png, file=paste("Network.HNPF.",NAME[i],".png",sep=""), height=800, width=3000)
dev.off()} else{
dev.copy(png, file=paste("Network.HNPF.",NAME[i],".png",sep=""), height=800, width=800)
dev.off()}
}
