## README folder structure of Cross-omics

# Data folder
* Results from the initial metabolite set creation and extension are saved in the 2019-09-12 folder. These are present here and not in the results, because they are used as existing data for the Gene_Metabolite_Mapping part of Cross-omics.
* All patient meta-data is present in the 'Crossomics_DBS....xlsx' and '.RData' files. The most recent date is used in the most recent cross-omics run.
* Patient metabolite Z-scores are present in the 'Project ...' folders.
* An excel file which includes all patients that have multiple disease genes is present as 'Patients_with....xlsx'.
* A list of all unique disease that are determined to be 'true' disease genes is present in the 'disease_genes.txt' file.
* A list of all genes used for the metabolite set creation and simulated WES results is present in 'All_Genes_....txt', with the version of Ensembl in the name.
* A list of all Gene-name to HMDB code translations is present in the 'HGNC2Uniprot....RDS' file
* Metabolite - reaction data from Recon3D is present.

# Results folder
* Results from different runs are present under the date they have been run. The most recent one (2019-12-10) is used for the paper as it is now (2020-01-09).
* The simulated WES results as used for the metabolite mapping are present in the Mock_genes folder. Inside the most recent date was used for the most recent run. 

# Plots folder
* All plots that I made since september. I used to have the R script create a folder by the date I made the plots, but I switched in December to have the date of the Cross-omics run as date of the folder. So all plots used in the paper are present in 2019-12-10.

# src_Metabolite_Set_Creation
All code used for creating and extending the gene-metabolite sets. 
* All scripts in the main folder are 'primary' scripts and are the ones that need to be run. 
* Scripts used by the main script are in the 'Supportive' folder.
* Scripts to run the R scripts from the terminal are present in the 'Run_files' folder. 

# src_Metabolite_Mapper
All code used for running the cross-omics scripts when the metabolite sets have been made and the patient Z-score files are present.
* All scripts in the main folder are 'primary' scripts and are the ones that need to be run. 
* Scripts used by the main script are in the 'Supportive' folder.
* Scripts to run the R scripts from the terminal are present in the 'Run_files' folder. 

# src_After_Analysis
All code used for plotting and summarizing results
* Concatenate_MSEA_results.R, followed by Summarise_MSEA_Results.R must be run to collate and summarise all the results coming from the metabolite mapper scripts before it can be visualised. These scripts will create .RDS output files in the Results folder
* Prioritised_gene_analysis.R is the script that I used for the majority of plotting (try-outs) and contains the code for the facet plot that visualises the fraction of correctly prioritised disease genes per parameter.
* Gene_met_set_stats and plotting.R summarise the metabolite set creation results and plot them.
* Nmr_aberrant_mets_per_patient.R summarises the number of aberrant metabolites for the disease gene of a patient per Z-score threshold. This script will create a .RDS output file in the Results folder.
* Distribution_Crossomics_vs_Random.R speaks for itself, but calculates and plots the fraction of 'correctly prioritizated' disease genes per parameter combination and what is expected by random.

# src_Other
Scripts that don't necessarily fall within one of the other categories, however, most of them must be run at some point.
* Between the creation of metabolite sets and running the metabolite mapper, CreateMockGeneSet.R, Patientxls_to_RData.R, Convert_xlsx_to_RDS.R, TranslateToHMDB.R and FixHMDB.R must be run.
1. CreateMockGeneSet.R speaks for itself. The number of genes in the Mock_gene sets and the number of sets must be adjust manually inside the script.
2. Patientxls_to_RData.R and Convert_xlsx_to_RDS.R are scripts make the excel files obsolete and convert them to more easily readable files for R. The former is used for transforming the Patient-meta-data and the latter to transform project patient Z-score files. In both scripts, the files to transform must be manually supplied.
3. TranslateToHMDB.R and FixHMDB.R must be run in that order to make the output of the metabolite set creation scripts more readily usable by converting as many kegg and chebi codes to HMDB and fix any of the resulting HMDB codes (I've had one particular metabolite always being translated to a faulty code which has a number too few).
* After running the MetaboliteMapper scripts, Missed_Patient_Seeds.R should be run to check which runs failed (on the HPC) and create files to rerun those.
* CheckDiseaseGenes.R is a script to check whether there are any disease genes present in the simulated WES results (Mock_genes), which shouldn't be the case, but if patient data has changed, you can check whether the mock-gene sets are still usable or not, which can be done with Exclude_seeds_with_disease_genes.R, but I don't recomment using this one, just remake the Mock_genes.
* CreatePatientSubset.R is not used currently, but could be used to create select only a subsection of patients to run.
* Zip_files.sh is just a simple example how to zip the results coming from the Cross-omics metabolite mapper run on the HPC. 

# Other_files
Older description files used to point out what all the code of the metabolite set creation and extention does in a way that I hope people understand.