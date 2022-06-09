library(dplyr)
library(readr)
library(tidyverse)
  # reading in, filtering, testing
full_GTDB_bac120 <- read.delim("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/bac120_metadata_r202.tsv", header = TRUE)
filt_GTDB_bac120 <- dplyr::filter(full_GTDB_bac120, mimag_high_quality == "t")
test_bac120 <- filt_GTDB_bac120$mimag_high_quality
full_GTDB_ar122 <- read.delim("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/ar122_metadata_r202.tsv", header = TRUE)
filt_GTDB_ar122 <- dplyr::filter(full_GTDB_ar122, mimag_high_quality == "t")
test_ar122 <- filt_GTDB_ar122$mimag_high_quality
complete <- bind_rows(filt_GTDB_bac120, filt_GTDB_ar122)
#print(test)
  # write TSV of filtered data
#write_tsv(filt_GTDB_bac120, file = "filtered_GTDB_bac120.tsv")
#write_tsv(filt_GTDB_ar122, file = "filtered_GTDB_ar120.tsv")
#write_tsv(complete, file = "filtered_GTDB_complete.tsv")
  # read in pre-filtered matrix

GTDB_genomes <- read.delim("/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/filtered_GTDB_genomes.txt", header = TRUE)

  # test presence or absence of fasta file paths
genomes_avail <- dplyr::filter(GTDB_genomes, GTDB_genomes$fasta_file_path != "NA")
