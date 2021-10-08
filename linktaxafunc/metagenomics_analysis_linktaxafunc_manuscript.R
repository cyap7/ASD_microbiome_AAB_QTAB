#==============================================================================
#
# Species-specific differential gene abundance - cleaning
#
#==============================================================================

# OVERVIEW:
# - for each MG_SampleID, iteratively read in files 
#	- for Romboutsia.timonensis should be n = 56+117 = 173
#		- subtract 1 as removing 4494902/BBC6589 = 172
# 	- for Erysipelatoclostridium.sp003024675, n=37
# - check if there is a Romboutsia column
# - if there is Romboutsia column, extract ID + Romboutsia column
# - rename the Romboutsia column to MG_SampleID

#==============================================================================
# Arguments
#==============================================================================

keep_id_dir 	<- commandArgs(trailingOnly = TRUE)[1]
data_dir 		<- commandArgs(trailingOnly = TRUE)[2]
suffix			<- commandArgs(trailingOnly = TRUE)[3]
extract			<- commandArgs(trailingOnly = TRUE)[4]
id_coln			<- commandArgs(trailingOnly = TRUE)[5]
diet_dir		<- commandArgs(trailingOnly = TRUE)[6]
func_count_dir	<- commandArgs(trailingOnly = TRUE)[7]
func_id_coln	<- commandArgs(trailingOnly = TRUE)[8]
out_dir			<- commandArgs(trailingOnly = TRUE)[9]
filen			<- commandArgs(trailingOnly = TRUE)[10]

#==============================================================================
# Libraries
#==============================================================================

library(data.table)
library(plyr)
library(tidyverse)
source("/90days/uqcyap3/software/ANCOM-v2.1/scripts/ancom_v2.1.R")

#==============================================================================
# Directories
#==============================================================================


# Directory and file names for the by_species data (containing 1 count table (species x genes) per sample)
# - file names of count tables within this directory must include the sample IDs (sample IDs will be grep-ed from the file list here)
data_dir <- "/PATH/TO/processed_data/autism_MCPv2_MGDBv2_functional/functional_profiles/by_species"
# - suffix to add to the file name to find and read in (included multiple options here depending on which gene level you want to look at)
suffix <- ".MICROBA.species.tsv.gz"
#suffix <- ".EC.species.tsv.gz"
#suffix <- ".MetaCyc_pathway.species.tsv.gz"
#suffix <- ".MetaCyc_group.species.tsv.gz"
# - colunm header corresponding to gene/feature of interest (included multiple options here depending on which gene level you want to look at)
id_coln <- "ID" # Microba genes, TCDB and EC
#id_coln <- "Pathway" # MetaCyc pathway
#id_coln <- "Group" # MetaCyc group

# Filters
extract <- "Romboutsia" # String with which to grep feature of interest. Genus "Romboutsia" is fine in this case as there are no other Romboutsia species in the dataset
keep_id_dir <- "/PATH/TO/aab_mgsample_Romboutsia.timonensis.id" # list of sample IDs who have Romboutsia timonensis present

# Pheno file
diet_dir <- "/PATH/TO/Diet/AES_ACRC_QTAB/AES_ACRC_QTAB_diet_PC_Microba_248_pcenergy_PC13_clr.csv" # Phenotype/dietary data file, which undergoes some cleaning in the metagenomics_analysis_linktaxafunc_ancom_direct.R and metagenomics_analysis_linktaxafunc_ancom_indirect.R scripts

# Directory and file name pointing to the count table of samples x genes
# - the directory
func_count_dir0 <- "/PATH/TO/processed_data/autism_MCPv2_MGDBv2_functional/functional_profiles/by_sample"
# - the file name (included multiple options here depending on which gene level you want to look at)
func_count_file <- "MICROBA.samples.tsv"
#func_count_file <- "EC.samples.tsv"
#func_count_file <- "MetaCyc_pathway.samples.tsv"
#func_count_file <- "MetaCyc_group.samples.tsv"
# - column name containing the gene/feature of interest
func_id_coln <- "VariableID"

# Output
# - directory for output
out_dir <- "/30days/uqcyap3/ASD/Data/3_metagenomics/diffabundance"
# - file name for output (included multiple options here depending on which gene level you want to look at)
filen <- "aab_Romboutsia.timonensis_MICROBAgenes"
#filen <- "aab_Romboutsia.timonensis_EC"
#filen <- "aab_Romboutsia.timonensis_MetaCyc_pathway"
#filen <- "aab_Romboutsia.timonensis_MetaCyc_group"


#==============================================================================
# Read in data
#==============================================================================

func_count_dir <- paste(func_count_dir0, func_count_file, sep = "/")

keep_id <- read.delim(keep_id_dir, header = F)
diet <- fread(diet_dir, header = T)
func_count <- fread(func_count_dir, header = T)

#==============================================================================
# Clean diet data
#==============================================================================

# Sort the pheno/covariates file
diet <- diet[order(diet$MG_SampleID),]

# Remove excluded individuals
diet <- diet %>% filter(!IID %in% c(4494902)) # Smith Magenis syndrome
diet <- diet %>% filter(participant_type %in% c("ASD", "SIB", "UNR")) %>% filter(participant_type_study %in% c("ASD", "SIB", "UNR", "UNR_QTAB")) # a double check, in case I eventually exclude others using a similar coding system

# Clean BSC data
diet$bristol_stool_chart_regroup <- NA
diet$bristol_stool_chart_regroup[which(diet$bristol_stool_chart %in% 1:2)] <- 2
diet$bristol_stool_chart_regroup[which(diet$bristol_stool_chart == 3)] <- 3
diet$bristol_stool_chart_regroup[which(diet$bristol_stool_chart == 4)] <- 4
diet$bristol_stool_chart_regroup[which(diet$bristol_stool_chart %in% 5:7)] <- 5
diet$bristol_stool_chart_regroup <- as.factor(diet$bristol_stool_chart_regroup)

# Make meds_antibiotic_current == 0 for QTAB (ie. just none reported for sensitivity analysis). Makes data munging easier (NAs are otherwise excluded)
diet$meds_antibiotic_current[which(is.na(diet$meds_antibiotic_current))] <- 0

# Get a diet file without missing data
# diet_full <- diet %>% filter(!is.na(PC1_diet_pe))

#==============================================================================
# Extract data
#==============================================================================

print("Extracting data")

#-----------------------------------
# First pass
#-----------------------------------

print("... First pass")

out <- NULL

# Do for the first version so not joining against empty data frame
i=1

print(keep_id[i,1])

# Read in file
file_ls <- list.files(data_dir)
in_file_n <- file_ls[grep(keep_id[i,1], file_ls)]
in_file_n <- in_file_n[grep(suffix, in_file_n)]
#in_file_n <- paste(data_dir, "/", keep_id[i,1], suffix, sep = "")
in_file <- fread(paste(data_dir, "/", in_file_n, sep = ""))

# Extract columns
check <- colnames(in_file)[grep(extract, colnames(in_file))]
ex <- c(id_coln, check)
tmp <- in_file[,..ex]

# Prune down file to non-zero gene abundances
tmp <- tmp[which(tmp[,2:ncol(tmp)] != 0)]

# Rename the column to be the keep_id name (MG_SampleID)
colnames(tmp) <- c("ID", keep_id[i,1])

# Join with out data frame
out <- tmp

#-----------------------------------
# All others
#-----------------------------------

print("... all others")

for (i in 2:nrow(keep_id)) {

	print(keep_id[i,1])
	# Read in file
	in_file_n <- file_ls[grep(keep_id[i,1], file_ls)]
	in_file_n <- in_file_n[grep(suffix, in_file_n)]
	in_file <- fread(paste(data_dir, "/", in_file_n, sep = ""))
	#in_file_n <- paste(data_dir, "/", keep_id[i,1], suffix, sep = "")
	#in_file <- fread(in_file_n)

	# Only proceed to extract if there is a extractable column!
	check <- colnames(in_file)[grep(extract, colnames(in_file))]
	if (length(check) != 0) {

		# Extract columns
		ex <- c(id_coln, check)
		tmp <- in_file[,..ex]

		# Prune down file to non-zero gene abundances
		tmp <- tmp[which(tmp[,2:ncol(tmp)] != 0)]

		# Rename the column to be the keep_id name (MG_SampleID)
		colnames(tmp) <- c("ID", keep_id[i,1])

		# Join with out data frame
		out <- join(out, tmp, by = "ID", type = "full")

	}
}

#==============================================================================
# Cleaning
#==============================================================================

print("Cleaning")
# Replace NA with 0
out <- as.matrix(out)
out[which((is.na(out)))] <- 0
out <- data.frame(out)
out <- data.frame(ID = out[,1], sapply(out[,2:ncol(out)], as.numeric))

#==============================================================================
# Outputs ("direct" analysis)
#==============================================================================

print("Write extracted tables")

# Write aggregated table
fwrite(out, paste(out_dir, "/", filen, ".tsv", sep = ""), 
	quote = F, sep = "\t")

# Write list of genes associated with the target
write.table(data.frame(out$ID), paste(out_dir, "/", filen, ".genes" , sep = ""), 
	quote = F, sep = "\t", row.names = F, col.names = F)

#==============================================================================
# Extract these genes from the overall set ("indirect" analysis)
#==============================================================================

print("Extract genes from the overall set")

# Table of genes that are expressed by Romboutsia --> look across all samples
func_count.df <- data.frame(func_count)
keep.tmp <- diet$MG_SampleID[which(diet$MG_SampleID %in% colnames(func_count.df))]
func_count.df <- data.frame(func_count.df[,1], func_count.df[,keep.tmp])
colnames(func_count.df)[1] <- "VariableID"
func_count2 <- func_count.df[which(func_count.df$VariableID %in% out$ID),]

# Write table
write.table(func_count2, paste(out_dir, "/", filen, "_genes_allMGIDs.tsv" , sep = ""), 
	quote = F, sep = "\t", row.names = F, col.names = T)

