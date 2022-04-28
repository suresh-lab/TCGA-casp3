# TCGA-casp3

Code and raw data for the genomic analyses for our pre-print found here: https://www.biorxiv.org/content/10.1101/2021.11.30.470656v1

This is the master branch which contains the original set of analysis, including the main analysis examining the effect of MK2 on survival in both NSCLC and other cancers. 

The "Revision 5" branch contains more extensive data, including: a) raw data (as rda or txt files) for the pan cancer analyses; b) additional pan cancer ananalyses performed during revisions of the original manuscript prior to publication; c) compressed count matrices for the various cancer types considered in our pan cancer analysis; d) count matrices and additional raw data used for the normal vs. tumor comparisons. That branch is stored in LFS format given the size of the repository, and will require lfs installation to git so that the files can be decompressed upon repo pull. 
General Notes:

scripts/
- main Cox PH and survival analyses, and pdfs of supplemental data
- scripts/old.versions: older versions of the main analysis md.
- MK2analysis_validation - this is the rmd and PDF file for the EA dataset

rawdata/
- includes raw data used in the scripts/ folder above. 
- Data Sources: OncoLnc, cBioPortal, TCGA (via TCGAbiolinks or via GDC portal)

TCGA-MK2_bootstrappedHR/
- contains scripts to generate the distribution of HR using all genes in LUAD
- APIaccess.rmd downloads a rda with the transcripts 
- In later analyses, we found it was easier just to grab this from firebrowse.org
- However, for this analysis, we did the preprocessing (preprocessing.R) locally
- LUAD_runbootstrap.R does most of the heavy lifting. Specifically, RunModelWithGene() is a function that runs a pre-specified Cox PH model by substituting in, sequentially, every gene in a large count matrix, thus giving a distribution of HRs. The results for each model, including model metrics, 
- LUAD_ProcessHR.R takes the results from LUAD_runbootstrap.R and plots the nice pretty figure in the paper. 

TCGA-MK2_panca/
- This was fun to write! So here, we go and get clinical data from a bunch of cancers(!).
- ./rawdata/MasterMK2data.csv contains the concatenated clinical data for a variety of tumor types
- ./scripts/PanCA_survival.analysis and annotation markdowns contain the code used to generate HRs for MK2  across cancer types
- So, the crux of the work is done by the PanCA_survival.analysis markdown returnModel() function
- We feed portions of the large concatenated dataset (one cancer at a time) to a function that calculates a Cox PH model on that subsetted data and returns the model results as a data structure, including metrics. We then extract from that list of lists the HR, CIs
- PanCA_survival.analysis rev2 contains HR (with CIs) but also the model metrics. 
