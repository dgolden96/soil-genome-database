# For GOLD data from 20 Jun, 2022
# Downloaded from https://gold.jgi.doe.gov/downloads

library(readxl)
library(tidyverse)
library(stringr)

# Zoey section: reads in data frame of JGI download data

gold_data_readme <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx", sheet = 1)
gold_data_study <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx", sheet = 2)
gold_data_biosample <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 3)
gold_data_organism <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 4)
gold_data_sequencing <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 5)
gold_data_analysis <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 6)

ecosystem_paths <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/GOLDs5levelEcosystemClassificationPaths.xlsx")

soil_paths <- ecosystem_paths %>% filter(grepl("Soil", `ECOSYSTEM TYPE`))
soil_path_ids <- soil_paths %>% select(`ECOSYSTEM PATH ID`) %>% unlist()

# Filter organisms to organisms found in soil
soil_organisms <- gold_data_organism %>% filter(`ORGANISM ECOSYSTEM PATH ID` %in% soil_path_ids)
# 16213 organisms
soil_organism_ids <- soil_organisms %>% select(`ORGANISM GOLD ID`) %>% unique() %>% unlist()

# Just testing the removal of fungi
no_fungi <- soil_organisms %>% filter(`ORGANISM NCBI SUPERKINGDOM` != "Eukaryota")

# Filter to genomes (analysis projects) from soil organisms
soil_sequencing <- gold_data_sequencing %>% filter(`ORGANISM GOLD ID` %in% soil_organism_ids)
soil_sequencing_ids <- soil_sequencing %>% select(`PROJECT GOLD ID`) %>% unlist()

# Filter to genomes (sequencing metadata) from soil organisms
# This file does not have all the completeness columns we want :(
soil_ap <- gold_data_analysis %>% filter(`AP PROJECT GOLD IDS` %in% soil_sequencing_ids)

genome_info <- merge(soil_ap, soil_sequencing, by.x = "AP PROJECT GOLD IDS", by.y = "PROJECT GOLD ID")
# Most of these columns can be removed, either before or after merging


# Next: add accession info as a column

accession_info <- genome_info$`AP GENBANK`

accession_info_test <- "[{\"genbankId\":\"ARFV00000000\",\"assemblyAccession\":\"GCA_000374445.1\"}]"
get_accession <- function(accession_info_test) {
	accession_info <- jsonlite::parse_json(accession_info_test)
	if (length(accession_info)==0) {
		return(NA)
	}
	out <- accession_info[[1]]
	print(out$assemblyAccession)
	to_return <- ifelse(is.null(out$assemblyAccession), NA, out$assemblyAccession)
	return(to_return)
}
get_accession(accession_info_test)

# Approach 1 - working
genome_info$NCBI_accession <- lapply(genome_info$`AP GENBANK`, get_accession)

# Approach 2 - working
genome_info$NCBI_accession <- NA
for (i in 1:nrow(genome_info)){
	print(i)
	genome_info$NCBI_accession[i] <- get_accession(genome_info$`AP GENBANK`[i])
}

# Check why ~1000 are NA
accession_NAs <- genome_info[is.na(genome_info$NCBI_accession),]

# Read in genbank and refseq data
gbrs_df_read <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/genbank_refseq_with_accession.csv")
gbrs_accession <- gbrs_df_read$X..assembly_accession

# Read in Dan's JGI query data for sanity check
JGI_dan_read <- read_tsv("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/JGI_data.tsv")
JGI_accession_dan <- JGI_dan_read$accession_col
zoey_dan_overlap <- intersect(genome_info$NCBI_accession, JGI_accession_dan)

# Filter downloaded JGI data to include only genomes for which we have gbrs accession numbers
# Also filter gbrs data to include only the genomes that are found in the JGI data
zoey_gbrs_overlap <- intersect(genome_info$NCBI_accession, gbrs_accession)
genome_info_gbrs_overlap <- dplyr::filter(genome_info, genome_info$NCBI_accession %in% zoey_gbrs_overlap)
gbrs_filt <- dplyr::filter(gbrs_df_read, gbrs_df_read$X..assembly_accession %in% genome_info_gbrs_overlap$NCBI_accession)
names(gbrs_filt)[names(gbrs_filt)=="X..assembly_accession"] <- "NCBI_accession"

# everything above this needs to be run at once before running the merges below
# sometimes you otherwise get an error I don't recognize

# Join selected JGI data with gbrs data on accession number, organism name, etc.
merge_genome_info <- merge(genome_info_gbrs_overlap, gbrs_filt, by = "NCBI_accession")
merge_genome_info <- merge(merge_genome_info, gold_data_organism, by.x = "AP NAME", by.y = "ORGANISM NAME")
merge_genome_info <- merge(merge_genome_info, gbrs_df_read, by.x = "NCBI_accession", by.y="X..assembly_accession")

# Select the fields we want and/or need
genome_info_select <- dplyr::select(merge_genome_info, NCBI_accession, 'AP NAME', 'AP PROJECT GOLD IDS', 'AP GOLD ID',
                                    'ORGANISM GOLD ID.x', 'ORGANISM NCBI TAX ID', taxid.x, species_taxid.x, organism_name.x,
                                    fasta_file_paths.x, 'ORGANISM NCBI GENUS', 'ORGANISM NCBI SPECIES')

# Change the names for clarity and Struo2 use
names(genome_info_select)[names(genome_info_select)=="species_taxid.x"] <- "ncbi_species_taxid"
names(genome_info_select)[names(genome_info_select)=="taxid.x"] <- "ncbi_taxid_nonspecies"
names(genome_info_select)[names(genome_info_select)=="fasta_file_paths.x"] <- "fasta_file_path_URL"
names(genome_info_select)[names(genome_info_select)=="NCBI_accession"] <- "partial_accession"
# Need to construct taxonomy column based on species column
for (i in 1:nrow(genome_info_select)){
  genome_info_select$'ORGANISM NCBI SPECIES'[i] <- sub(" ", "; ", genome_info_select$'ORGANISM NCBI SPECIES'[i])
}
names(genome_info_select)[names(genome_info_select)=="ORGANISM NCBI SPECIES"] <- "ncbi_taxonomy_partial"
fasta_URLs <- genome_info_select$fasta_file_path_URL
genome_names_JGI <- substr(fasta_URLs, start=83, stop=121)
genome_names_JGI <- paste0(genome_names_JGI, "_genomic.fna.gz")
genome_names_JGI_trunc <- str_extract(string = fasta_URLs, pattern = "(?<=v1/).*(?=.fna.gz)")
genome_names_JGI_fixed <- paste0(genome_names_JGI_trunc, ".fna.gz")
working_fasta <- c()

# this will take ~2 hours to run. It downloads the .fna.gz files for all the genomes in JGI meeting our criteria
  #for (i in 1:length(fasta_URLs)) {
  #  if (genome_names_trunc_fixed[i] == "NA.fna.gz") next
  #  tryCatch({
  #    download.file(fasta_URLs[i], paste0("/projectnb2/microbiome/dgolden/Struo2/data/JGI_genomes/", genome_names_trunc_fixed[i]))
  #    append(working_fasta, paste0("/projectnb2/microbiome/dgolden/Struo2/data/JGI_genomes/", genome_names_trunc_fixed[i]))
  #  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #}

  # let's extend accession from "GCA_#####" to "GB_GCA_#####"
JGI_accession <- genome_info_select$partial_accession
for (i in 1:length(JGI_accession)) {
  if (substr(JGI_accession[i], start=1, stop=3) == "GCA") {
    JGI_accession[i] <- paste0("GB_", JGI_accession[i])
  }
  if (substr(JGI_accession[i], start=1, stop=3) == "GCF") {
    JGI_accession[i] <- paste0("RS_", JGI_accession[i])
  }
}
genome_info_select$full_accession <- JGI_accession

taxonomy_species <- genome_info_select$ncbi_taxonomy_partial
taxonomy_species <- gsub(";", "", taxonomy_species)
taxonomy_species <- paste0("START", taxonomy_species)
taxonomy_species <- paste0(taxonomy_species, "END")
taxonomy_species <- str_extract(string = taxonomy_species, pattern = "(?<=START).*(?=END)")
taxonomy_species <- paste0("s__", taxonomy_species)
taxonomy_genus <- genome_info_select$ncbi_taxonomy_partial
taxonomy_genus <- paste0("START", taxonomy_genus)
taxonomy_genus <- str_extract(string = taxonomy_genus, pattern = "(?<=START).*(?=; )")
taxonomy_genus <- paste0("g__", taxonomy_genus)
full_tax <- paste(taxonomy_genus, taxonomy_species, sep = ";")
genome_info_select$ncbi_taxonomy_full <- full_tax

JGI_fasta_filepaths <- paste0("/projectnb2/microbiome/dgolden/Struo2/data/JGI_genomes/", genome_names_JGI_fixed)
genome_info_select$fasta_file_path <- JGI_fasta_filepaths
names(genome_info_select)[names(genome_info_select)=="AP NAME"] <- "ncbi_organism_name"
genome_info_select <- dplyr::filter(genome_info_select, genome_info_select$fasta_file_path != "/projectnb2/microbiome/dgolden/Struo2/data/JGI_genomes/NA.fna.gz")
organism_names_JGI <- genome_info_select$ncbi_organism_name
unique_organism_names <- make.unique(organism_names_JGI)
genome_info_select$unique_ncbi_organism_name <- unique_organism_names
#write_tsv(genome_info_select, "/projectnb2/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_corrected.tsv")

genome_info_distinct_name <- distinct(genome_info_select, ncbi_organism_name, .keep_all = TRUE)
genome_info_distinct_name_acc <- distinct(genome_info_distinct_name, full_accession, .keep_all = TRUE)
genome_info_distinct_accession <- distinct(genome_info_select, full_accession, .keep_all = TRUE)
genome_info_distinct_id <- distinct(genome_info_select, ncbi_species_taxid, .keep_all = TRUE)

unique_names_and_accession <- make.unique(genome_info_distinct_accession$ncbi_organism_name)
genome_info_distinct_accession$unique_ncbi_organism_name <- unique_names_and_accession

#write_tsv(genome_info_distinct_accession, "/projectnb2/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_distinct_accession_unique_names.tsv")
#write_tsv(genome_info_distinct_name, "/projectnb2/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_distinct_name.tsv")

premade_tsv <- read_tsv("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release202/taxdump/taxdump/taxID_info.tsv")
premade_accession <- premade_tsv$name
JGI_accession_distinct_name <- genome_info_distinct_name_acc$full_accession
JGI_premade_overlap <- intersect(JGI_accession_distinct_name, premade_accession)
JGI_accession_filtered <- JGI_accession_distinct_name[! JGI_accession_distinct_name %in% JGI_premade_overlap]

genome_info_distinct_accession_premade_filter <- dplyr::filter(genome_info_distinct_accession, full_accession %in% JGI_accession_filtered)
#write_tsv(genome_info_distinct_accession, "/projectnb2/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_distinct_name_acc_premade_filter.tsv")

