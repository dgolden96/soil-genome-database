# For GOLD data from 20 Jun, 2022
# Downloaded from https://gold.jgi.doe.gov/downloads

library(readxl)
library(tidyverse)

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

gbrs_df_read <- read.csv("/projectnb2/microbiome/dgolden/Struo2/genomes_for_file_paths/genbank_refseq_with_accession.csv")
gbrs_accession <- gbrs_df_read$X..assembly_accession
JGI_final_read <- read_tsv("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/JGI_data.tsv")
JGI_accession <- JGI_final_read$accession_col
zoey_dan_overlap <- intersect(genome_info$NCBI_accession, JGI_accession)
zoey_gbrs_overlap <- intersect(genome_info$NCBI_accession, gbrs_accession)
genome_info_gbrs_overlap <- dplyr::filter(genome_info, genome_info$NCBI_accession %in% zoey_gbrs_overlap)
merge_genome_info <- merge(genome_info_gbrs_overlap, gold_data_organism, by.x = "AP NAME", by.y = "ORGANISM NAME")
merge_genome_info <- merge(merge_genome_info, gbrs_df_read, by.x = "NCBI_accession", by.y="X..assembly_accession")
genome_info_select <- dplyr::select(merge_genome_info, NCBI_accession, 'AP NAME', 'AP PROJECT GOLD IDS', 'AP GOLD ID',
                                    'ORGANISM GOLD ID.x', bioproject, biosample, taxid, species_taxid, organism_name,
                                    fasta_file_paths, 'ORGANISM NCBI GENUS', 'ORGANISM NCBI SPECIES')
names(genome_info_select)[names(genome_info_select)=="species_taxid"] <- "taxID_col"
  # most likely NCBI species taxid
names(genome_info_select)[names(genome_info_select)=="taxID_col"] <- "ncbi_taxid"
names(genome_info_select)[names(genome_info_select)=="fasta_file_paths"] <- "fasta_file_path_col"
names(genome_info_select)[names(genome_info_select)=="NCBI_accession"] <- "accession"
# Need to construct taxonomy column based on species column
for (i in 1:nrow(genome_info_select)){
  genome_info_select$'ORGANISM NCBI SPECIES'[i] <- sub(" ", "; ", genome_info_select$'ORGANISM NCBI SPECIES'[i])
}
names(genome_info_select)[names(genome_info_select)=="ORGANISM NCBI SPECIES"] <- "ncbi_taxonomy"
#write_tsv(genome_info_select, "/projectnb2/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_corrected.tsv")
