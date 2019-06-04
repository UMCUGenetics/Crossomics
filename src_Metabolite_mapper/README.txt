README - Crossomics Metabolite Mapper

02-04-2019 - Marten Kerkhofs



This repository is for all Crossomics Metabolite Mapper scripts. The main script (for now) is GeneMetaboliteMapper_Marten_AdductSums.R, which is copied and modified from the similarly named GeneMetaboliteMapper_..._.R files. 

Files directly dependent on the main file are getPValues_Marten.R and performMSEA_Marten.R. The getPValues file must be renamed in the near future to represent the fact that it only calculates Z-scores, no P.values.

The other files have unknown use for now.

It is made to run on my Mac and some document-paths are hardcoded inside it. 
Before being able to run it at all, the Metab drive must be loaded.