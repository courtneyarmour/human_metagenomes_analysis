# human_metagenomes_analysis
A Metagenomic Meta-Analysis Reveals Functional Signatures of Health and Disease in the Human Gut Microbiome

File descriptions:
# make_circos_karyotype.R
Code for creating the karyotype used to make the circos plot in the manuscript

# make_dataframes_for_cplm.R
Used to produce a dataframe for each kegg object in each disease with the KEGG object abundance and sample disease status. 
Used for running regression models.

# manuscript_figures.R
Code for creating the figures in the manuscript (except the circos plot)

# metagenome_analysis_functions.R
Functions created for this project, this file is sourced in most of the other scripts. 

# ms_review.R
This file contains the code for the analyses described in Supplemental text 1 of the manuscript. 

# run_model.R
The code for running CPGLM on the dataframes created with the make_dataframes_for_cplm.R script.
