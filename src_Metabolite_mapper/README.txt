README - Crossomics Metabolite Mapper

2019-09-04 - Marten Kerkhofs


IN GENERAL:

This repository is for all necessary Crossomics Metabolite Mapper scripts. The data used and the mock genes are in the parent folders (Crossomics/Data and Crossomics/Results, the latter for the mock genes). Intensity values for the different metabolomics projects can be found on the Metab (Y:) drive, but I copy them to my Crossomics/Data.


HOW TO USE:

To run the MetaboliteMapper on the HPC: All files necessary for running these scripts need to be transferred to the HPC, with the file structures:
1. src_Metabolite_Mapper/... files:
	a. GeneMetaboliteMapper.R and its supportive scripts (see PIPELINE HIERARCHY/STRUCTURE)
	b. run_seeds.sh, run_array.sh and run_MetabMapper.sh. These files will call each other (in this order). Any parameter that the script uses can be changed in one of these files (see .SH FILES).
2. Data/... files:
	a. xxx.RData file containing all patients, their disease genes, data-file locations and patient numbers
	b. Adductsums files (Project xxx/whatever/path/follows/Bioinformatics/adductsums files): Patient/control metabolomics intensity files, which can be copied from the Metab (Y:) drive.
	c. gene specific metabolite sets. There are many in a folderstructure like this: yyyy-mm-dd/maxrxnxx/mss_x_HMDBtranslated/gene_whatever.RDS, where x needs to be replaced by a valid parameter condition.
	d. All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt: a file containing all genes downloaded from Ensembl with 4 columns: Geen stable ID, Gene type, Gene name and HGNC ID;
	e. mets2remove.RDS: a file containing all metabolites that should be removed before MSEA analysis. It contains 5 columns: Name, chebi, kegg, pubchem and inchi. Basically all of them are metabolite names in different formats.
3. Results/Mock_genes/mockgenefiles.txt: simple text files containing a certain number of mock genes. The file name tells the number of mock genes and the seed used for sampling them. see OTHER SCRIPTS for more info.


.SH FILES:
** run_seeds.sh runs run_array.sh and gives it 1 seed. All used seeds need to be changed/entered inside the run_seeds script and these MUST be the same as those used for creating the mock_gene files (its in their name).
** run_array.sh will call on run_MetabMapper.sh and submit that job via qsub. The directory in which all scripts are present must be correctly present in this script as 'code_dir'. Any changes in reserving HPC space and time must be adjusted in this script.
** run_MetabMapper.sh will run the GeneMetabMapper.R script via R. There are 2 important things with this script: 
	1. The location of R must be present in this script (I use R 3.6.0, while the standard HPC R is different). If running the scripts malfunctions, it is a good idea to also check whether all R libraries are correctly installed on the HPC. 
	2. All parameters must be adjusted/entered in this script as it will loop over all of them. maxrxn and step are the literal parameter values, but threshold is only the NUMBER of parameter options (as this parameter actually consists of 2 values inside the .R scripts and this was easier to code).


IMPORTANT NOTES:

All excel files used for these scripts must be converted to .RData and .RDS files. This is to ensure the scripts can work without any java-dependencies and make it easier to be HPC-compatible.


PIPELINE HIERARCHY/STRUCTURE:

The main script is GeneMetaboliteMapper.R. Scripts used by the main script are in the Supportive folder and are: getZValues_11042019.R, MSEA.R, sourceDir.R and Generate_av_Z_scores.R. The excel script 'genExcelFileShort.R' has been removed in one of the recent versions and replaced by making .RData files.
Other scripts do not have any use in the Metabolite_Mapper script and are moved to 'old/unnecessary' folders.


OTHER SCRIPTS

** The CreateMockGeneSet.R script creates mock gene text files with the seed and number of mock genes in the file title.
** CreatePatientSubset.R is not used atm, but works the same way as CreateMockGeneSet, but then for patients
** TranslateToHMDB.R is a script that simply translates all chebi and kegg codes within the metabolite set files to HMDB. It should be run locally after creating the metabolite sets to ensure that the metabolite mapper is able to be run without Java.


