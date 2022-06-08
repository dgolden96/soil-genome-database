library(dplyr)
library(readr)
  # reading in, filtering, testing
full_GTDB <- read.delim("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/bac120_metadata_r202.tsv", header = TRUE)
filt_GTDB <- dplyr::filter(full_GTDB, full_GTDB$mimag_high_quality == "t")
test <- filt_GTDB$mimag_high_quality
#print(test)
  # write TSV of filtered data
#write_tsv(filt_GTDB, file = "filtered_GTDB.tsv")

  # read in pre-filtered matrix
GTDB_genomes <- read.delim("/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/filtered_GTDB_genomes.txt", header = TRUE)

  # test presence or absence of fasta file paths
genomes_avail <- dplyr::filter(GTDB_genomes, GTDB_genomes$fasta_file_path != "NA")
