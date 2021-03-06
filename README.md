# README: ASD_microbiome_AAB_QTAB code
#### Author: Chloe X Yap
#### Date: 8 October 2021

[![DOI](https://zenodo.org/badge/376450308.svg)](https://zenodo.org/badge/latestdoi/376450308)

## Description:
This repository contains code used in the AAB and QTAB stool metagenomics analysis, and accompanying the paper titled, **"Autism-gut microbiome connections may reflect autism-related dietary diversity"**.

## Contents:
`metagenomics_analysis_taxonomic_manuscript.Rmd`

* contains code for most analyses (eg. microbiome diversity, differential abundance, linear models)
* N.B. the code contained here is specific to the taxonomic for demonstration purposes.
	However, where mentioned in the manuscript, identical analyses were performed for functional datasets.

`diet/`
* `AES_PCs_ACRC_QTAB_manuscript.html`
	* .html knitted .Rmd document
	* includes code and figures from exploratory analysis of dietary data

`osca_variancecomponent/`
* `metagenomics_analysis_osca_prep_manuscript.Rmd`
	* .Rmd file used to pre-process microbiome count data (taxonomic and functional).
	* output of this code is used as input into OSCA, using commands in metagenomics_analysis_osca_manuscript.sh

* `metagenomics_analysis_osca_manuscript.sh`
	* .sh file with OSCA commands to generate ORMs, and to perform variance component analysis

`linktaxafunc/`
* `metagenomics_analysis_linktaxafunc_manuscript.R`
	* R script to extract:
		* count table for genes directly and specifically mapping to Romboutsia timonensis ("direct" analysis)
	 	* count table for genes that are within the Romboutsia timonensis genome, but can be encoded by any other species ("indirect" analysis)
* `metagenomics_analysis_linktaxafunc_ancom_direct_manuscript.R`
	* R code used to perform the "direct" analysis
	* input the "direct" analysis output from `metagenomics_analysis_linktaxafunc_manuscript.R` (eg. "aab_Romboutsia.timonensis_MICROBAgenes.txt")
	* differential abundance analysis (ANCOMv2.1)

##### Updated: 14 May 2022

Updated `metagenomics_analysis_osca_prep_manuscript.Rmd`
- to perform the clr-transform, we used the mixOmics::logratio.transfo() function. We inadvertently had used the transpose of the intended matrix. We have re-run the OREML analyses and confirmed that the impact on the presented results are minior, and that the conclusions stand as in the original manuscript. We provide here the updated version of Fig.2 which directly pertains to these results. We sincerely apologise for any inconvenience caused by this error.