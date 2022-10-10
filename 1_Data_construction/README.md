# Note to generate ASV tables
Run F1_dada2_ITS_IonTorrent.r and P1_dada2_16S_IonTorrent.r. 
Before starting, please install the fastq files from XXXX, and also download database (silva_nr_v138_train_set.fa.gz and sh_general_release_dynamic_all_10.05.2021.fasta), and change the lines with "# CHANGE ME" to your own directry.

Because our IonTorrent setting changed between the runs to generate the fastq files for this project, some of the fastq files had full forward primer sequence but other files lost the first one or two bases from the forward primer sequence. Therefore, these R scripts contains the process to remove frist one or two bases from the samples with full forward primer sequences to equalize the start base of the forward primer sequences in all the samples.
In bacterial 16S RNA fastq files, the fastq files namesd as from "...16S_001_R..." to "...16S_010_R...", from "...16S_012_R..." to "...16S_025_R...", from "...16S_029_R..." to "...16S_036_R...", and from "...16S_043_R..." to "...16S_049_R..." have the full forward primer sequence, but the other files does not have the first two bases GT.
In fungal ITS fastq files, the fastq files namesd as from "...ITS_103_R..." to "...ITS_114_R..." have the full forward primer sequence, but the other files does not have the first one bases T.

# Data availability
As mentioned above, the sequence data is available online. Please contact me if you need any other input files (e.g., metadata.csv).

# File order
For ITS, please run F1 -> F2 -> F3 files. For 16S rRNA, please run P1 -> P2 -> ... -> P4 files.

