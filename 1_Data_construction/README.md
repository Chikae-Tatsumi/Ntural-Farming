# To obtain the ASV files
Run dada2_ITS_IonTorrent.r and dada2_16S_IonTorrent.r.
Before starting, please install the fastq files from XXXX, and also download database (silva_nr_v138_train_set.fa.gz and sh_general_release_dynamic_all_10.05.2021.fasta), and change the lines with "# CHANGE ME" to your own directry.

Because our IonTorrent setting changed between the runs to generate the fastq files for this project, some of the fastq files had full forward primer sequence but other files lost the first one or two bases from the forward primer sequence. Therefore, in these R scripts, we removed frist one or two bases from the samples with full forward primer sequences to equalize the start base of the forward primer sequences in all the samples.
