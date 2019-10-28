README - Crossomics Metabolite Mapper

2019-09-04 - Marten Kerkhofs


IN GENERAL:

This repository is for all necessary Crossomics Metabolite Mapper scripts. The data used and the mock genes are in the parent folders (Crossomics/Data and Crossomics/Results, the latter for the mock genes). Intensity values for the different metabolomics projects can be found on the Metab (Y:) drive, but I copy them to my Crossomics/Data.

Important note on creating mock gene sets (see OTHER SCRIPTS)


HOW TO USE:

To run the MetaboliteMapper on the HPC: All files necessary for running these scripts need to be transferred to the HPC, with the file structures:
1. src_Metabolite_Mapper/... files:
	a. GeneMetaboliteMapper.R and its supportive scripts (see PIPELINE HIERARCHY/STRUCTURE)
	b. run_full_MetabMapper.sh and run_R.sh. These files will call each other (in this order). Any parameter that the script uses can be changed in one of these files (see .SH FILES).
2. Data/... files:
	a. xxx.RData file containing all patients, their disease genes, data-file locations and patient numbers
	b. Adductsums files (Project xxx/whatever/path/follows/Bioinformatics/adductsums files): Patient/control metabolomics intensity files, which can be copied from the Metab (Y:) drive.
	c. gene specific metabolite sets. There are many in a folderstructure like this: yyyy-mm-dd/maxrxnxx/mss_x_HMDBtranslated/gene_whatever.RDS, where x needs to be replaced by a valid parameter condition.
	d. All_Genes_Ensembl_apr_2019_GRCh38p12_extended.txt: a file containing all genes downloaded from Ensembl with 4 columns: Geen stable ID, Gene type, Gene name and HGNC ID;
	e. mets2remove.RDS: a file containing all metabolites that should be removed before MSEA analysis. It contains 5 columns: Name, chebi, kegg, pubchem and inchi. Basically all of them are metabolite names in different formats.
3. Results/Mock_genes/mockgenefiles.txt: simple text files containing a certain number of mock genes. The file name tells the number of mock genes and the seed used for sampling them. see OTHER SCRIPTS for more info.


.SH FILES:
** run_full_MetabMapper.sh creates a qsub job of run_R.sh and gives it 1 file name containing a seed, all parameter combinations with each one being 1 single string and some directories where to find my own R version and the coding directory. All used seeds (their filenames) are read in from a folder with this kind path "/hpc/shared/dbg_mz/marten/Crossomics_2019_10_25/Results/Mock_genes". It is vital that the different parameter values are separated with a ',' (comma) inside the string, so my R script can easily separate them. It is also vital to supply the number of patients that must be run; this is because the job-submission works on an array-basis, which takes simple 1-n sub-jobs and gives that number to the script which is called.

** run_R.sh will run the GeneMetabMapper.R script via Rscript. Everything inside this script is given to it via run_full_MetabMapper.sh, so the only thing that could possibly be changed here is when the R script has had its name changed or the number or order or parameter values which are supplied to it has changed.


IMPORTANT NOTES:

All excel files used for these scripts must be converted to .RData and .RDS files. This is to ensure the scripts can work without any java-dependencies and make it easier to be HPC-compatible.


PIPELINE HIERARCHY/STRUCTURE:

The main script is GeneMetaboliteMapper.R. Scripts used by the main script are in the Supportive folder and are: getZValues_11042019.R, MSEA.R, sourceDir.R and Generate_av_Z_scores.R. The excel script 'genExcelFileShort.R' has been removed in one of the recent versions and replaced by making .RData files.
Other scripts do not have any use in the Metabolite_Mapper script and are moved to 'old/unnecessary' folders.


OTHER SCRIPTS

** The CreateMockGeneSet.R script creates mock gene text files with the seed and number of mock genes in the file title. IMPORTANT NOTE: these need to be remade every time a new set of patients is used. The patient-data is one of the inputs and this script makes sure that none of the disease genes is present in the mock gene sets.
** CreatePatientSubset.R is not used atm, but works the same way as CreateMockGeneSet, but then for patients
** TranslateToHMDB.R is a script that simply translates all chebi and kegg codes within the metabolite set files to HMDB. It should be run locally after creating the metabolite sets to ensure that the metabolite mapper is able to be run without Java.


