# To obtain the ASV files
Run dada2_ITS_IonTorrent.r and dada2_16S_IonTorrent.r.
Before starting, please install the fastq files from XXXX, and also download database (silva_nr_v138_train_set.fa.gz and sh_general_release_dynamic_all_10.05.2021.fasta), and change the lines with "# CHANGE ME" to your own directry.

Because our IonTorrent setting changed between the runs to generate the fastq files for this project, some of the fastq files had full primer sequence but other sequences lost the first one or two bases. Therefore, in these R scripts, we removed frist one or two bases from the samples with full primer sequences to equalize the start base in all the samples.
