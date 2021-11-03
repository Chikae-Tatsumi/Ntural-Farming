# DADA2 Pipeline Tutorial (1.12)
# See http://benjjneb.github.io/dada2/tutorial.html

# Getting ready
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

DATABASE = "~/R/Database/silva_nr_v138_train_set.fa.gz" # CHANGE ME to the directory containing the database

# For with GT
setwd("~/R/Analysis/4_Paddy/HNPF/16S/withGT")  ## CHANGE ME to the directory containing the fastq files.
filez <- list.files()
file.rename(from=filez, to=sub(pattern=".fastq", replacement=".fastq.gz", filez))

fnFs <- sort(list.files(getwd(), pattern = ".fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(fnFs[1:2])

# Filter and trim
filtFs_withGT <- file.path(getwd(), "filtN", paste0(sample.names, ".fastq.gz"))
names(filtFs_withGT) <- sample.names

out_withGT <- filterAndTrim(fnFs, filtFs_withGT, maxN = 0, maxEE = 2, truncLen = 240,
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = 19)  # on mac, set multithread = TRUE
head(out_withGT)



# For without GT
setwd("~/R/Analysis/4_Paddy/HNPF/16S/withoutGT")  ## CHANGE ME to the directory containing the fastq files.
filez <- list.files()
file.rename(from=filez, to=sub(pattern=".fastq", replacement=".fastq.gz", filez))

# Define which is forward fastq and reverse fastq
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(getwd(), pattern = ".fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(fnFs[1:2])

# Filter and trim
filtFs_withoutGT <- file.path(getwd(), "filtN", paste0(sample.names, ".fastq.gz"))
names(filtFs_withoutGT) <- sample.names

out_withoutGT <- filterAndTrim(fnFs, filtFs_withoutGT, maxN = 0, maxEE = 2, truncLen = 238,
    truncQ = 2, minLen = 48, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = 17)  # on mac, set multithread = TRUE
head(out_withoutGT)


# Merging with GT and without GT (Chikae's original)
filtFs <- c(filtFs_withGT, filtFs_withoutGT)
out <- rbind(out_withGT, out_withoutGT)

out <- out[order(rownames(out)),]
filtFs <- filtFs[order(names(filtFs))]

setwd("~/R/Analysis/4_Paddy/HNPF/16S")  ## CHANGE ME to the directory containing the fastq files.

# If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails
# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) #plot

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab) 
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #>90% is good?

# Track reads through the pipeline
# look at the number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,file="track.txt")
write.table(seqtab.nochim, file="seqtabnochim.txt")

# Assign taxonomy
# Install "silva_nr_v132_train_set.fa.gz" from: https://zenodo.org/record/1172783#.XUmvQ_ZFw2w
# Other taxonomic reference database: http://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim,DATABASE, multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
# If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) and see if this fixes the assignments. 

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file="taxonomy.txt")

samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxa))
#store the DNA sequences of our ASVs in the refseq slot of the phyloseq object

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Remove Mitochondria, Chloroplast, and Eukaryote
ps_removed = subset_taxa(ps,(
                             Family  != "Mitochondria"|is.na(Family) &
                             Class   != "Chloroplast"|is.na(Class)  &
                             Kingdom  != "Eukaryota" ))

# To output ASV table
otu_table.t<-t(ps_removed@otu_table)
ps.t<-cbind(otu_table.t,ps_removed@tax_table)
write.table(ps.t,  file="ASV_table.txt")

# Rarefication
ps.rarefied = rarefy_even_depth(ps_removed, rngseed=1, sample.size=min(sample_sums(ps_removed)), replace=F)
otu_table.t<-t(ps.rarefied@otu_table)
ps.t<-cbind(otu_table.t,ps.rarefied@tax_table)
write.table(ps.t,  file="rarefied_ASV_table.txt")
