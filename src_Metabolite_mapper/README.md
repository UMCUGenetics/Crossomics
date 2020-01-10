# README - Crossomics Metabolite Mapper

2019-09-04 - Marten Kerkhofs
2020-01-10 - Marten Kerkhofs


## IN GENERAL:
This repository is for all necessary Crossomics Metabolite Mapper scripts. The data used and the mock genes are in the parent folders (Crossomics/Data and Crossomics/Results, the latter for the mock genes). Intensity values for the different metabolomics projects can be found on the Metab (Y:) drive, but I copied them to my Crossomics/Data folder.

Important note on creating mock gene sets (see OTHER SCRIPTS)


## HOW TO USE:
To run the MetaboliteMapper on the HPC: All files necessary for running these scripts need to be transferred to the HPC, with the file structures:
1. src_Metabolite_Mapper/... files:
	a. GeneMetaboliteMapper.R and its supportive scripts in a 'Supportive' folder (see PIPELINE HIERARCHY/STRUCTURE)
	b. run_full_MetabMapper.sh and run_R.sh. These files will call each other (in this order). Any parameter that the script uses can be changed in one of these files (see Run_files).
2. Data/... files:
	a. xxx.RData file containing all patients, their disease genes, data-file locations and patient numbers. The naming looks like this at the moment: Crossomics_DBS_Marten_trimmed....RData and then select the most recent one.
	b. Project-specific Z-score files (Project xxx/whatever/path/follows/Bioinformatics/projectxxx.RDS files): Patient/control metabolomics intensity files, which can be copied from the Metab (Y:) drive. In the Metab drive, they are present as .xlsx files. However, for faster run times, I convert them to .RDS files, which can simply be done via the script: Convert_xlsx_to_RDS.R (in the "../src_Other folder")
	c. gene specific metabolite sets. There are many in a folderstructure like this: yyyy-mm-dd/maxrxnxx/mss_x_HMDBtranslated/gene_whatever.RDS, where x needs to be replaced by a valid parameter condition.
	d. All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt: a file containing all genes downloaded from Ensembl with 4 columns: Geen stable ID, Gene type, Gene name and HGNC ID.
	e. mets2remove.RDS: a file containing all metabolites that should be removed before MSEA analysis. It contains 5 columns: Name, chebi, kegg, pubchem and inchi. Basically all of them are metabolite names in different formats.
3. Results/Mock_genes/\[mockgenefiles\].txt: simple text files containing a certain number of mock genes. The file name tells the number of mock genes and the seed used for sampling them. see OTHER SCRIPTS for more info.

## Starting and testing the script
* To actually start/run the whole thing, if everything is in place, you simply need to start it via ```./run_full_MetabMapper.sh```. 
* This can also be tested beforehand with the help of these two scripts: 1: `test_working_run_full_MetabMapper.sh` and 2: `test_time_run_full_MetabMapper.sh`. 
1. The former tests whether the scripts are all in place and everything works in general (it lets the main script perform just a few iterations to run through the whole thing). 
2. The latter can be used to give an indication for the time and memory that is necessary for running the whole thing on the HPC. To make an educated guess as to what space and time is necessary, the whole run can be made faster by adjusting those requirements inside the run_full_MetabMapper.sh. To do this, the testing script performes a few full runs and the output is emailed to the user (the email adress is hardcoded in).


## Run_files:
Beforehand: the .sh scripts can email when it they are ready or error. This emailadress is still mine (Marten's) so you need to change that in the .sh calls to qsub.
* run_full_MetabMapper.sh creates a qsub job of run_R.sh and gives it 1 file name containing a seed, all parameter combinations in a single string and some directories where to find my own R version and the coding directory. All used seeds (their filenames) are read in from a folder with this kind path "/hpc/shared/dbg_mz/marten/Crossomics_2019_10_25/Results/Mock_genes" (but now I removed everything from the HPC and it is present in /hpc/dbg_mz/users/Marten/Crossomics...). It is vital that the different parameter values are separated with a ',' (comma) inside the string, so my R script can easily separate them. It is also vital to supply the number of patients that must be run; this is because the job-submission works on an array-basis, which takes simple 1-to-n sub-jobs and gives that number to the script which is called.
* run_R.sh will run the GeneMetabMapper.R script via Rscript. Everything inside this script is given to it via run_full_MetabMapper.sh, so the only thing that could possibly be changed here is when the R script has had its name changed or the number or order or parameter values which are supplied to it has changed.


## IMPORTANT NOTES:
All excel files used for these scripts must be converted to .RData and .RDS files. This is to ensure the scripts can work without any java-dependencies, making it easier to be HPC-compatible and to improve the speed of the scripts.


## PIPELINE HIERARCHY/STRUCTURE:
The main script is GeneMetaboliteMapper.R. Scripts used by the main script are in the Supportive folder and are: MSEA.R, get_Z_score_matrix.R and sourceDir.R.


## OTHER SCRIPTS (present in ../src_Other)
These scripts can be run (most need to be run at least once) after the creation of metabolite sets and before running the metabolite mapper.
* The CreateMockGeneSet.R script creates mock gene text files with the seed and number of mock genes in the file title. IMPORTANT NOTE: these need to be remade every time a new set of patients is used. The patient-data is one of the inputs and this script makes sure that none of the disease genes is present in the mock gene sets.
* CreatePatientSubset.R is not used atm, but works the same way as CreateMockGeneSet, but then for patients
* TranslateToHMDB.R is a script that simply translates all chebi and kegg codes within the metabolite set files to HMDB. It should be run locally after creating the metabolite sets to ensure that the metabolite mapper is able to be run without Java on the HPC.
* Convert_xlsx_to_RDS.R obviously converts xlsx to RDS files. It is specifically used to convert the Z-score / intensity files that are made by the general metabolomics pipeline. It is a very basic script, but it needs you to change the working directory every time you use it to convert whatever file you want from, and to, the same folder as the working directory.

