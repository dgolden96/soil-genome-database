library(tidyverse)
library(dplyr)
library(readr)

genbank_assembly <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/assembly_summary_genbank.csv", header = TRUE)
mod_vec <- genbank_assembly$ftp_path
genbank_assembly$fasta_file_paths <- mod_vec
genbank_assembly$fasta_file_paths <- paste0(genbank_assembly$fasta_file_paths, "_genomic.fna.gz")
refseq_assembly <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/assembly_summary_refseq.csv", header = TRUE)
mod_vecRS <- refseq_assembly$ftp_path
refseq_assembly$fasta_file_paths <- mod_vecRS
refseq_assembly$fasta_file_paths <- paste0(refseq_assembly$fasta_file_paths, "_genomic.fna.gz")
#write_tsv(genbank_assembly, "genbank_with_filepaths.tsv")
#write_tsv(refseq_assembly, "refseq_with_filepaths.tsv")
gb_accession <- genbank_assembly$X..assembly_accession
rs_accession <- refseq_assembly$X..assembly_accession
gbrs_accession <- c(gb_accession, rs_accession)
gbrs_df <- rbind(genbank_assembly, refseq_assembly)
full_tsv <- read_tsv("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/metadata/filtered_GTDB_complete.tsv")
full_accession <- full_tsv$accession
full_accession_trim <- substr(full_accession, start=4, stop=18)
overlap <- intersect(gbrs_accession, full_accession_trim)
full_tsv$accession_trim <- full_accession_trim
filtered_GTDB_with_FASTA <- dplyr::filter(full_tsv, full_tsv$accession_trim %in% overlap)
#write_tsv(filtered_GTDB_with_FASTA, "filt_GTDB_FASTA.tsv")

JGI_org_df <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/MAG_organisms.csv")
JGI_proj_df <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/MAG_projects.csv")
org_ids_names <- JGI_org_df$Organism.Organism.Name
proj_ids_names <- JGI_proj_df$Analysis.Project.Analysis.Project.Name
names(JGI_org_df)[names(JGI_org_df)=="Organism.Organism.Name"] <- "samples_col"
names(JGI_proj_df)[names(JGI_proj_df)=="Analysis.Project.Analysis.Project.Name"] <- "samples_col"
#print(org_ids_names == proj_ids_names)
  # this confirms that the organism and project names are the same
JGI_org_proj <- merge(JGI_org_df, JGI_proj_df, by="samples_col")
  # join on same names, similar to SQL joins
#write_tsv(JGI_org_proj, "JGI_data.tsv")
JGI_accession <- JGI_org_proj$Analysis.Project.Assembly.Accession
JGI_overlap <- intersect(gbrs_accession, JGI_accession)
JGI_overlap_df <- dplyr::filter(JGI_org_proj, Analysis.Project.Assembly.Accession %in% JGI_overlap)
gbrs_JGI_overlap_df <- dplyr::filter(gbrs_df, X..assembly_accession %in% JGI_overlap)
names(JGI_org_proj)[names(JGI_org_proj)=="Analysis.Project.Assembly.Accession"] <- "accession_col"
names(gbrs_JGI_overlap_df)[names(gbrs_JGI_overlap_df)=="X..assembly_accession"] <- "accession_col"
JGI_final <- merge(JGI_org_proj, gbrs_JGI_overlap_df, by="accession_col")
names(JGI_final)[names(JGI_final)=="Organism.NCBI.Taxonomy.ID"] <- "taxID_col"
names(JGI_final)[names(JGI_final)=="fasta_file_paths"] <- "fasta_file_path_col"
