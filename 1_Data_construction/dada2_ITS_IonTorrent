# DADA2 ITS Pipeline Workflow (1.8)
# See https://benjjneb.github.io/dada2/ITS_workflow.html

# Getting Ready
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

DATABASE = "~/R/Database/sh_general_release_dynamic_all_10.05.2021.fasta" # CHANGE ME to the directory containing the database

# Identify primers
FWD <-	"CCGTAGGTGAACCTGCGG" ## ITS1 "TCCGTAGGTGAACCTGCGG"
REV <-	"GCTGCGTTCTTCATCGATGC" ## ITS2

setwd("~/R/Analysis/4_Paddy/HNPF/ITS/withT")  ## CHANGE ME to the directory containing the fastq files.

fnFs <- sort(list.files(getwd(), pattern = ".fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)
filtFs <- file.path("~/R/Analysis/4_Paddy/HNPF/ITS_database=all", paste0(sample.names, ".fastq"))
out <- filterAndTrim(fnFs, filtFs, trimLeft = 1)  # on mac, set multithread = TRUE

setwd("~/R/Analysis/4_Paddy/HNPF/ITS")  ## CHANGE ME to the directory containing the fastq files.
filez <- list.files()
file.rename(from=filez, to=sub(pattern=".fastq", replacement=".fastq.gz", filez))

fnFs <- sort(list.files(getwd(), pattern = ".fastq", full.names = TRUE))

# to ensure we have the right primers, and the correct orientation of the primers on the reads
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name ))
head(sample.names)
fnFs.filtN <- file.path(getwd(), "filtN", paste0(sample.names, ".fastq.gz"))

#to “pre-filter” the sequences just to remove those with  ambiguous bases (Ns)
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

# Identifying and counting the primers on one set of paired end FASTQ files(Sufficient, because all the files were created using the same library preparation)
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]) )

#Remove Primers
cutadapt <- "/miniconda2/bin/cutadapt" 
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(getwd(), "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnFs)) # dummy file

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2,
                            "-o", fnFs.cut[i],  # output file                             
                            fnFs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = ".fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2,
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

# see https://github.com/benjjneb/dada2/issues/384
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")

rownames(track) <- sample.names
head(track)
write.table(track,file="track.txt")

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, DATABASE, multithread = TRUE, tryRC = TRUE)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file="taxonomy.txt")
rownames(seqtab.nochim) <- sample.names
write.table(seqtab.nochim, file="seqtabnochim.txt")

samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxa))
# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Remove others than fungi
ps_removed = subset_taxa(ps,(Kingdom  == "k__Fungi" ))
                             
# To output OTU table
otu_table.t<-t(ps_removed@otu_table)
ps.t<-cbind(otu_table.t,ps_removed@tax_table)
write.table(ps.t,  file="ASV_table.txt")

# Rarefication
ps.rarefied = rarefy_even_depth(ps_removed, rngseed=1, sample.size=min(sample_sums(ps_removed)), replace=F)
otu_table.t<-t(ps.rarefied@otu_table)
ps.t<-cbind(otu_table.t,ps.rarefied@tax_table)
write.table(ps.t,  file="rarefied_ASV_table.txt")
