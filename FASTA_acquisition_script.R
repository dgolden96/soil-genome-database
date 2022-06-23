library(tidyverse)
library(dplyr)
library(readr)

genbank_assembly <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/assembly_summary_genbank.csv", header = TRUE)
mod_vec <- genbank_assembly$ftp_path
genome_names <- substr(mod_vec, start=58, stop=95)
genome_names <- paste0(genome_names, "_genomic.fna.gz")
genome_names <- paste0("/", genome_names)
full_paths <- paste0(mod_vec, genome_names)
genbank_assembly$fasta_file_paths <- full_paths
refseq_assembly <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/assembly_summary_refseq.csv", header = TRUE)
mod_vecRS <- refseq_assembly$ftp_path
genome_namesRS <- substr(mod_vecRS, start=58, stop=95)
genome_namesRS <- paste0(genome_namesRS, "_genomic.fna.gz")
genome_namesRS <- paste0("/", genome_namesRS)
full_pathsRS <- paste0(mod_vecRS, genome_namesRS)
refseq_assembly$fasta_file_paths <- full_pathsRS

#write.csv(genbank_assembly, "/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/genbank_with_filepaths.csv")
#write.csv(refseq_assembly, "/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/refseq_with_filepaths.csv")
gb_accession <- genbank_assembly$X..assembly_accession
rs_accession <- refseq_assembly$X..assembly_accession
gbrs_accession <- c(gb_accession, rs_accession)
gbrs_df <- rbind(genbank_assembly, refseq_assembly)
full_tsv <- read_tsv("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/filtered_GTDB_complete.tsv")
full_accession <- full_tsv$accession
full_accession_trim <- substr(full_accession, start=4, stop=18)
overlap <- intersect(gbrs_accession, full_accession_trim)
full_tsv$accession_trim <- full_accession_trim
filtered_GTDB <- dplyr::filter(full_tsv, full_tsv$accession_trim %in% overlap)
filtered_GTDB_with_FASTA <- merge(filtered_GTDB, gbrs_df, by.x = "accession_trim", by.y = "X..assembly_accession")
filtered_GTDB_with_FASTA_select <- dplyr::select(filtered_GTDB_with_FASTA, accession, accession_trim, checkm_completeness, checkm_contamination,
                                                 gtdb_genome_representative, gtdb_taxonomy, mimag_high_quality, ncbi_bioproject,
                                                 ncbi_genbank_assembly_accession, ncbi_organism_name, ncbi_species_taxid, ncbi_taxid,
                                                 ncbi_taxonomy, bioproject, biosample, taxid, species_taxid, organism_name, assembly_level,
                                                 ftp_path, fasta_file_paths)
#write_tsv(filtered_GTDB_with_FASTA, "/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/filt_GTDB_FASTA.tsv")
#write_tsv(filtered_GTDB_with_FASTA_select, "/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/filt_GTDB_FASTA_select.tsv")
#write_csv(gbrs_df, "/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/genbank_refseq_with_accession.csv")
