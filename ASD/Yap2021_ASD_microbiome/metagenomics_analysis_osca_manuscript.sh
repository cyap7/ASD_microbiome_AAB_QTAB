#==============================================================================
#
# Metagenomics OSCA analysis (taxonomic families with median count > 0)
#
#==============================================================================

# LATER: mgrm analysis with rare taxa too

data=/PATH/TO/DATA/DIRECTORY
out=/PATH/TO/OUTPUT/DIRECTORY

# Taxonomic datasets file prefixes
type=taxa
taxlev=species_common_clrcheck
taxlev=genus_common_clrcheck
taxlev=family_common_clrcheck
taxlev=species_rm0_clrcheck
taxlev=species_rm0_rare_clrcheck
taxlev=species_rm0_rare_01
taxlev=metaphlan2species_common_clrcheck

# Functional dataset file prefixes
type=func
taxlev=EClevel4_common_clrcheck
taxlev=EClevel4_rm0_rare_clrcheck
taxlev=EClevel4_rm0_rare_01
taxlev=EClevel3_common_clrcheck
taxlev=TCDB_common_clrcheck
taxlev=MetaCycpathway_common_clrcheck
taxlev=Microba_common_clrcheck
taxlev=Microba_rm0_rare_clrcheck
taxlev=Microba_rm0_rare_01
taxlev=humann2_common_clrcheck

# Dietary dataset file prefix
filen=AES_ACRC_QTAB_diet_food_245

#==============================================================================
#
# Generate the ORMs
#
#==============================================================================

#------------------------------------------------------------------------------
# Microbiome datasets
#------------------------------------------------------------------------------

# Prepare to generate the ORM the ORM (in OSCA)
# - for taxonomic and function microbiome datasets
~/osca \
--efile ${data}/${type}_count_${taxlev}_t.tsv \
--make-bod \
--out ${data}/${type}_count_${taxlev}_t 

# Generate ORM + PCA
~/osca \
--befile ${data}/${type}_count_${taxlev}_t  \
--make-orm-bin \
--orm-alg 2 \
--out ${data}/${type}_count_${taxlev}_t 

# Generate PCs of the ORM
~/osca \
--orm ${data}/${type}_count_${taxlev}_t  \
--pca 20 \
--out ${data}/${type}_count_${taxlev}_t 

#------------------------------------------------------------------------------
# Dietary matrix dataset
#------------------------------------------------------------------------------

# - for the dietary dataset (food frequency matrix with ordinal variables)
~/osca \
--efile ${data}/${filen}.tsv \
--make-bod \
--out ${data}/${filen}

# Generate ORM + PCA
~/osca \
--befile ${data}/${filen} \
--make-orm-bin \
--orm-alg 2 \
--out ${data}/${filen}

# Generate PCs of the ORM
~/osca \
--orm ${data}/${filen}  \
--pca 20 \
--out ${data}/${filen}

#==============================================================================
#
# Variance analysis: microbiome datasets (taxonomic + functional)
#
#==============================================================================

# The commands below is structured by phenotype
# ie. all ASD-related OREML analyses are presented first (primary + sensitivity),
# then dietary PC, dietary diversity, bristol stool chart OREML analyses, etc.

#----------------------------
# OREML: ASD as phenotype
#----------------------------
# Primary analysis
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_asd_nobristol_reml

# - sensitivity: no covariates
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--out ${out}/${type}_count_${taxlev}_t_asd_nocov_reml

# - sensitivity: no Abx
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--remove ${data}/antibiotics_current.id \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_asd_nobristol_noabx_reml

# - sensitivity: no QTAB, all ages <= 10
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--keep ${data}/ageu10_noqtab.id \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_asd_nobristol_ageu10noqtab_reml

# - sensitivity: no QTAB, all ages <= 10, no covariates
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--keep ${data}/ageu10_noqtab.id \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--out ${out}/${type}_count_${taxlev}_t_asd_nobristol_ageu10noqtab_nocov_reml

# - sensitivity: no SIB
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--keep ${data}/nosib.id \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--out ${out}/${type}_count_${taxlev}_t_asd_nobristol_nosib_reml

# - sensitivity: no SIB, no covariates
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--keep ${data}/nosib.id \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--out ${out}/${type}_count_${taxlev}_t_asd_nobristol_nosib_nocov_reml

#----------------------------
# Dietary PC as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_dietPC1_nogroup_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_dietPC2_nogroup_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_dietPC3_nogroup_reml

# Dietary PC (with ASD group as covariate)
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${type}_count_${taxlev}_t_dietPC1_inclgroup_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${type}_count_${taxlev}_t_dietPC2_inclgroup_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${type}_count_${taxlev}_t_dietPC3_inclgroup_reml

# - try pe_PC1 without age (as associated with age in linear regression)
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${type}_count_${taxlev}_t_dietPC1_inclgroup_noage_reml

#----------------------------
# Dietary diversity
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_foodshannon.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${type}_count_${taxlev}_t_foodshannon_inclgroup_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_foodshannon.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_foodshannon_nogroup_reml

#----------------------------
# Bristol Stool Chart as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_bristol_incldiet_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_bristol_nodiet_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--remove ${data}/antibiotics_current.id \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_bristol_nodiet_noabx_reml

#----------------------------
# Age as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--reml-maxit 1000 \
--out ${out}/${type}_count_${taxlev}_t_age_nodiet_reml

# - sensitivity: include diet
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_age.qcov \
--qcovar ${data}/osca_metagenomics_dietPC.qcov \
--reml-maxit 1000 \
--out ${out}/${type}_count_${taxlev}_t_age_incldiet_reml

# - sensitivity: no Abx
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--remove ${data}/antibiotics_current.id \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_age_nodiet_noabx_reml

# - sensitivity: no QTAB, age <= 10
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--keep ${data}/ageu10_noqtab.id \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_age_nodiet_ageu10noqtab_reml

# - sensitivity, no SIB
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--keep ${data}/nosib.id \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_age_nodiet_nosib_reml

# - this matrix is not invertible
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-maxit 1000 \
--out ${out}/${type}_count_${taxlev}_t_age_incldiet_inclgroup_reml

# BMI as phenotype
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_bmi.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_bmi_nodiet_reml

# - sensitivity: including dietary PCs
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_bmi.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_bmi_incldiet_reml

#----------------------------
# Sleep as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_sleep_cshq_total.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_sleep_nodiet_reml

# - sensitiivty: including diet
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_sleep_cshq_total.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_sleep_incldiet_reml

#----------------------------
# Development quotient as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_devq_wisc_msel_nih.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_devq_nodiet_reml

# - sensitivity: including diet
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_devq_wisc_msel_nih.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_devq_incldiet_reml

#----------------------------
# CD4 T cell (choose this as seem to be impt https://www.nature.com/articles/s41422-020-0332-7)
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_CD4Tcell.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${type}_count_${taxlev}_t_CD4Tcell_nobristol_reml

#==============================================================================
#
# Variance analyses: dietary matrix dataset(from AES food frequency score)
#
#==============================================================================

filen=AES_ACRC_QTAB_diet_food_245

#----------------------------
# ASD as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${filen}_t_asd_nobristol_nodiet_reml

#----------------------------
# Dietary PC including group as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${filen}_t_dietPC1_inclgroup_reml

~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${filen}_t_dietPC2_inclgroup_reml

~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--out ${out}/${filen}_t_dietPC3_inclgroup_reml

#----------------------------
# Dietary diversity (no diet covariates)
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_foodshannon.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${filen}_t_foodshannon_nogroup_reml

#----------------------------
# Bristol (no diet covariates)
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${filen}_t_bristol_nodiet_reml

#----------------------------
# Age as phenotype (no diet covariates)
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${filen}_t_age_nodiet_reml

#----------------------------
# BMI as phenotype (no diet covariates)
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_bmi.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${filen}_t_bmi_nodiet_reml

#----------------------------
# Sleep as phenotype (no diet covariates)
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_sleep_cshq_total.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${filen}_t_sleep_nodiet_reml

#----------------------------
# Development quotient as phenotype (no diet covariates)
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_devq_wisc_msel_nih.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--out ${out}/${filen}_t_devq_nodiet_reml

#==============================================================================
#
# Variance analyses without covariates: microbiome datasets
#
#==============================================================================

type=taxa
taxlev=species_common_clrcheck
taxlev=species_rm0_rare_clrcheck
taxlev=species_rm0_rare_01

type=func
taxlev=EClevel4_common_clrcheck
taxlev=TCDB_common_clrcheck
taxlev=MetaCycpathway_common_clrcheck
taxlev=Microba_common_clrcheck
taxlev=Microba_rm0_rare_clrcheck
taxlev=Microba_rm0_rare_01

#----------------------------
# ASD
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--out ${out}/${type}_count_${taxlev}_t_asd_nocov_reml

#----------------------------
# Dietary PC1-3
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--out ${out}/${type}_count_${taxlev}_t_dietPC1_nocov_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--out ${out}/${type}_count_${taxlev}_t_dietPC2_nocov_reml

~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--out ${out}/${type}_count_${taxlev}_t_dietPC3_nocov_reml

#----------------------------
# Dietary diversity
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_foodshannon.pheno \
--out ${out}/${type}_count_${taxlev}_t_foodshannon_nocov_reml

#----------------------------
# Bristol Stool Chart as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--out ${out}/${type}_count_${taxlev}_t_bristol_nocov_reml

#----------------------------
# Age as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_age.qcov \
--reml-maxit 1000 \
--out ${out}/${type}_count_${taxlev}_t_age_nocov_reml

#----------------------------
# BMI as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_bmi.pheno \
--out ${out}/${type}_count_${taxlev}_t_bmi_nocov_reml

#----------------------------
# Sleep as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_sleep_cshq_total.pheno \
--out ${out}/${type}_count_${taxlev}_t_sleep_nocov_reml

#----------------------------
# Development quotient as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${type}_count_${taxlev}_t \
--pheno ${data}/osca_metagenomics_devq_wisc_msel_nih.pheno \
--out ${out}/${type}_count_${taxlev}_t_devq_nocov_reml

#==============================================================================
#
# Variance analyses: dietary matrix dataset, no covariates
#
#==============================================================================

filen=AES_ACRC_QTAB_diet_food_245

#----------------------------
# ASD as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--out ${out}/${filen}_t_asd_nocov_reml

#----------------------------
# Dietary PC including group
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--out ${out}/${filen}_t_dietPC1_nocov_reml

~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--out ${out}/${filen}_t_dietPC2_nocov_reml

~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--out ${out}/${filen}_t_dietPC3_nocov_reml

#----------------------------
# Dietary diversity
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_foodshannon.pheno \
--out ${out}/${filen}_t_foodshannon_nocov_reml

#----------------------------
# Bristol
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--out ${out}/${filen}_t_bristol_nocov_reml

#----------------------------
# Age as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_age.qcov \
--out ${out}/${filen}_t_age_nocov_reml

#----------------------------
# BMI as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_bmi.pheno \
--out ${out}/${filen}_t_bmi_nocov_reml

#----------------------------
# Sleep as phenotype
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_sleep_cshq_total.pheno \
--out ${out}/${filen}_t_sleep_nocov_reml

#----------------------------
# Development quotient as phenotypE
#----------------------------
~/osca \
--reml \
--orm ${data}/${filen} \
--pheno ${data}/osca_metagenomics_devq_wisc_msel_nih.pheno \
--out ${out}/${filen}_t_devq_nocov_reml

#==============================================================================
#
# Multiple ORM analyses (ie. fit multiple ORMs in the model simultaneously)
#
#==============================================================================

# The hashed-out text below the vi commands list the matrices to include 
# within that .txt file.
# In turn, that .txt file is used as input into the multiple ORM analyses
# - eg. the first example
#	```vi ${data}/taxa_count_species_common_rare_clrcheck.txt```
# 	includes the names of 2 ORMs generated in the OSCA software package
#		1. common species
#		2. rare species

# Prepare mgrm files
vi ${data}/taxa_count_species_common_rare_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_species_rm0_rare_clrcheck_t

vi ${data}/taxa_count_species_common_rare_01_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_species_rm0_rare_01_t

vi ${data}/taxa_count_species_genus_family_common_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_genus_common_clrcheck_t
# /PATH/TO/taxa_count_family_common_clrcheck_t

vi ${data}/taxa_count_species_genus_family_common_species_rare_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_species_rm0_rare_clrcheck_t
# /PATH/TO/taxa_count_genus_common_clrcheck_t
# /PATH/TO/taxa_count_family_common_clrcheck_t

vi ${data}/taxa_count_species_common_rare_func_count_EClevel4_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_species_rm0_rare_clrcheck_t
# /PATH/TO/func_count_EClevel4_common_clrcheck_t

vi ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_species_rm0_rare_clrcheck_t
# /PATH/TO/func_count_EClevel4_common_clrcheck_t
# /PATH/TO/func_count_TCDB_common_clrcheck_t

vi ${data}/func_count_EClevel4_TCDB_common_clrcheck.txt
# /PATH/TO/func_count_EClevel4_common_clrcheck_t
# /PATH/TO/func_count_TCDB_common_clrcheck_t

vi ${data}/func_count_EClevel4_MetaCycpathway_common_clrcheck.txt
# /PATH/TO/func_count_EClevel4_common_clrcheck_t
# /PATH/TO/func_count_MetaCycpathway_common_clrcheck_t

vi ${data}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt
# /PATH/TO/func_count_EClevel4_common_clrcheck_t
# /PATH/TO/func_count_TCDB_common_clrcheck_t
# /PATH/TO/func_count_MetaCycpathway_common_clrcheck_t

vi ${data}/taxa_count_species_common_rare_func_count_EClevel4_MetaCycpathway_common_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_species_rm0_rare_clrcheck_t
# /PATH/TO/func_count_EClevel4_common_clrcheck_t
# /PATH/TO/func_count_MetaCycpathway_common_clrcheck_t

vi ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt
# /PATH/TO/taxa_count_species_common_clrcheck_t
# /PATH/TO/taxa_count_species_rm0_rare_clrcheck_t
# /PATH/TO/func_count_EClevel4_common_clrcheck_t
# /PATH/TO/func_count_TCDB_common_clrcheck_t
# /PATH/TO/func_count_MetaCycpathway_common_clrcheck_t

# Dietary PC1 and species (common + rare)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_clrcheck_t_dietPC1_mgrm_reml

# Dietary PC1 and species (common + rare with 0/1 encoding)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_01_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_01_clrcheck_t_dietPC1_mgrm_reml

# Add all collapsed taxonomic data (species, genus, family)
# - does this assume independence? Or no?
# - leave MetaCyc pathway out of dietary PC1 analyses as doesn't explain any variance
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_clrcheck_t_dietPC1_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_species_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_species_rare_clrcheck_t_dietPC1_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_common_clrcheck_t_dietPC1_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/func_count_EClevel4_TCDB_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC1.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/func_count_EClevel4_TCDB_common_clrcheck_t_dietPC1_mgrm_reml

# Dietary PC2 and species (common + rare)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_clrcheck_t_dietPC2_mgrm_reml

# Add all collapsed taxonomic data (species, genus, family)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_clrcheck_t_dietPC2_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_species_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC2.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_species_rare_clrcheck_t_dietPC2_mgrm_reml

# - left out MetaCyc as non-invertible matrix
# - species (common + rare) and EClevel4 + TCDB is non-vertible
# ~/osca \
# --reml \
# --multi-orm ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_clrcheck.txt \
# --pheno ${data}/osca_metagenomics_dietPC2.pheno \
# --qcovar ${data}/osca_metagenomics_age.qcov \
# --covar ${data}/osca_metagenomics_sex_group.cov \
# --reml-no-lrt \
# --reml-no-constrain \
# --out ${out}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_common_clrcheck_t_dietPC2_mgrm_reml

# - EC + TCDB +/- MetaCyc pathway analyses also had non-invertible matrix error
# ~/osca \
# --reml \
# --multi-orm ${data}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt \
# --pheno ${data}/osca_metagenomics_dietPC2.pheno \
# --qcovar ${data}/osca_metagenomics_age.qcov \
# --covar ${data}/osca_metagenomics_sex_group.cov \
# --reml-no-lrt \
# --reml-no-constrain \
# --out ${out}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck_t_dietPC2_mgrm_reml

# Dietary PC3 and species (common + rare)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_clrcheck_t_dietPC3_mgrm_reml

# Add all collapsed taxonomic data (species, genus, family)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_clrcheck_t_dietPC3_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_species_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_species_rare_clrcheck_t_dietPC3_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck_t_dietPC3_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_dietPC3.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--reml-no-constrain \
--out ${out}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck_t_dietPC3_mgrm_reml

# Bristol Stool Scale and species (common + rare)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_clrcheck_t_bristol_mgrm_reml

# Add all collapsed taxonomic data (species, genus, family)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_clrcheck_t_bristol_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_species_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_species_rare_clrcheck_t_bristol_mgrm_reml

# EClevel4 + MetaCyc
~/osca \
--reml \
--multi-orm ${data}/func_count_EClevel4_MetaCycpathway_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/func_count_EClevel4_MetaCycpathway_common_clrcheck_t_bristol_mgrm_reml

# Species (common + rare) + EC + TCDB + MetaCyc
# - failed to converge after 3000 iterations
# ~/osca \
# --reml \
# --multi-orm ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt \
# --pheno ${data}/osca_metagenomics_bristol.pheno \
# --qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
# --covar ${data}/osca_metagenomics_sex_group.cov \
# --reml-no-lrt \
# --reml-maxit 3000 \
# --out ${out}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck_t_bristol_mgrm_reml

# Species (common + rare) + EC + TCDB
# - failed to converge after 3000 iterations
# ~/osca \
# --reml \
# --multi-orm ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_clrcheck.txt \
# --pheno ${data}/osca_metagenomics_bristol.pheno \
# --qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
# --covar ${data}/osca_metagenomics_sex_group.cov \
# --reml-no-lrt \
# --reml-maxit 3000 \
# --out ${out}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_common_clrcheck_t_bristol_mgrm_reml

# Species (common + rare) + EC + MetaCyc
# - failed to converge after 3000 iterations
# ~/osca \
# --reml \
# --multi-orm ${data}/taxa_count_species_common_rare_func_count_EClevel4_MetaCycpathway_common_clrcheck.txt \
# --pheno ${data}/osca_metagenomics_bristol.pheno \
# --qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
# --covar ${data}/osca_metagenomics_sex_group.cov \
# --reml-no-lrt \
# --reml-maxit 3000 \
# --out ${out}/taxa_count_species_common_rare_func_count_EClevel4_MetaCycpathway_common_clrcheck_t_bristol_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_bristol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck_t_bristol_mgrm_reml

# Age and species (common + rare, no dietary PCs as covariates)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_clrcheck_t_age_mgrm_reml

# Add all collapsed taxonomic data (species, genus, family)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_clrcheck_t_age_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_genus_family_common_species_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_genus_family_common_species_rare_clrcheck_t_age_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck_t_age_mgrm_reml

~/osca \
--reml \
--multi-orm ${data}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck.txt \
--pheno ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex_group.cov \
--reml-no-lrt \
--reml-no-constrain \
--out ${out}/func_count_EClevel4_TCDB_MetaCycpathway_common_clrcheck_t_age_mgrm_reml

# BMI and species (common + rare, no dietary PCs as covariates)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_bmi.pheno \
--qcovar ${data}/osca_metagenomics_age.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_clrcheck_t_bmi_mgrm_reml

# ASD diagnosis and species (common + rare)
~/osca \
--reml \
--multi-orm ${data}/taxa_count_species_common_rare_clrcheck.txt \
--pheno ${data}/osca_metagenomics_casecontrol.pheno \
--qcovar ${data}/osca_metagenomics_age_dietPC.qcov \
--covar ${data}/osca_metagenomics_sex.cov \
--reml-no-lrt \
--out ${out}/taxa_count_species_common_rare_clrcheck_t_asd_mgrm_reml

