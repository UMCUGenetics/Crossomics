README - Crossomics Metabolite Mapper

2019-09-04 - Marten Kerkhofs



This repository is for all necessary Crossomics Metabolite Mapper scripts. The data used is in the parent folder (Data) and on the Metab drive in a project specific path. 

For an HPC version of the Metabolite Mapper, the data present in the Metab drive is copied to the Data folder. In this version, all uses of excel are also replaced with RData and RDS files to make the script Java-independent.

The main script is GeneMetaboliteMapper.R. Scripts used by the main script are in the Supportive folder and are: getZValues_11042019.R, MSEA.R, sourceDir.R and Generate_av_Z_scores.R. The excel script 'genExcelFileShort.R' has been removed in one of the recent versions and replaced by making .RData files.
Other scripts do not have any use in the Metabolite_Mapper script and are moved to 'old/unnecessary' folders.



