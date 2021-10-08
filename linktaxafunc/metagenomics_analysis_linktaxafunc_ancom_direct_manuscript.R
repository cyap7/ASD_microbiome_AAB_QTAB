#==============================================================================
#
# Species-specific differential gene abundance - ANCOM
#
#==============================================================================

#==============================================================================
# Arguments
#==============================================================================

out_dir			<- commandArgs(trailingOnly = TRUE)[1] # output directory
filen			<- commandArgs(trailingOnly = TRUE)[2] # eg. "aab_Romboutsia.timonensis_MICROBAgenes"
diet_dir		<- commandArgs(trailingOnly = TRUE)[3] # eg. "/PATH/TO/Diet/AES_ACRC_QTAB/AES_ACRC_QTAB_diet_PC_Microba_248_pcenergy_PC13_clr.csv"
filter_comm_dir	<- commandArgs(trailingOnly = TRUE)[4] # file including "common" species to include. This just helps to filter the count file and reduce memory
filter_rare_dir	<- commandArgs(trailingOnly = TRUE)[5] # file including "common" species to include. This just helps to filter the count file and reduce memory

#==============================================================================
# Libraries
#==============================================================================

library(data.table)
library(plyr)
library(tidyverse)
source("/PATH/TO/software/ANCOM-v2.1/scripts/ancom_v2.1.R")

#==============================================================================
# Read in data
#==============================================================================

out <- read.delim(paste(out_dir, "/", filen, ".tsv" , sep = ""), header = T)
diet <- fread(diet_dir, header = T)

if (filter_comm_dir != "NA" | filter_rare_dir != "NA") {
	filter_comm <- fread(filter_comm_dir)
	filter_rare <- fread(filter_rare_dir)
	
	out <- data.frame(out)
	colnames(out)[1] <- "VariableID"
	
	# Filter down into those with >10 non-zero values
	filter <- rbind(filter_comm, filter_rare)
	out <- out[which(out$VariableID %in% filter$V2),]
}

colnames(out)[1] <- "VariableID"

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
# Test for differences in genes among those with the target
#==============================================================================

print("ANCOM analysis: differences in genes directly from this target")

#------------------------------------------------------------
# Prepare data
#------------------------------------------------------------

print("... Prepare data for ANCOM")

# Test for differential abundance of genes, in individuals with this taxa
diet_extract <- diet[which(diet$MG_SampleID %in% colnames(out)),]
identical(diet_extract$MG_SampleID, colnames(out)[2:ncol(out)])

# Prepare data
metadata <- diet_extract %>% dplyr::select(c("MG_SampleID", "ASD", "proforma_sex", "age", "PC1_diet_pe", "PC2_diet_pe", "PC3_diet_pe", "MG_FID"))

metadata$ASD <- as.factor(metadata$ASD)
metadata$MG_FID[which(metadata$MG_FID == 0)] <- NA

ancom_taxa <- out[,c(2:ncol(out))]
ancom_taxa <- ancom_taxa[,which(colnames(ancom_taxa) %in% metadata$MG_SampleID)]
rownames(ancom_taxa) <- out$VariableID

# Check variables in the correct order
identical(colnames(ancom_taxa), diet_extract$MG_SampleID)

# Get dimensions of input data
dim(ancom_taxa)

ancom_table <- feature_table_pre_process(ancom_taxa, metadata, "MG_SampleID", group_var = NULL, out_cut = 0.05, zero_cut = 1, lib_cut = 0, neg_lb)

saveRDS(ancom_table, paste(out_dir, "/", filen, "_ancom_table.rds", sep = ""))

#------------------------------------------------------------
# ANCOM + cov: sex + age + dietary PCs
#------------------------------------------------------------

print("... Case-control analysis: covariates sex + age + dietary PC1-3")
# Run tests with covariates (age + sex + diet)
ancom_out_sex_age_diet <- ANCOM(ancom_table$feature_table, ancom_table$meta_data, ancom_table$structure_zeros, main_var = "ASD", p_adj_method = "BH", adj_formula = "proforma_sex + age + PC1_diet_pe + PC2_diet_pe + PC3_diet_pe")

#ancom_out_sex_age_diet$out <- ancom_out_sex_age_diet$out[order(ancom_out_sex_age_diet$out$W, decreasing = T),]

write.table(ancom_out_sex_age_diet$out, paste(out_dir, "/", filen, "_ANCOM_directgenes_sexagediet.txt", sep = ""),
	col.names = T, row.names = F, sep = "\t", quote = F)

saveRDS(ancom_out_sex_age_diet$fig, paste(out_dir, "/", filen, "_ANCOM_directgenes_sexagediet.rds", sep = ""))

