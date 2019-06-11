README - Crossomics Metabolite Set Creation

31-05-2019 - Marten Kerkhofs
Update: 07-06-2019 - Marten Kerkhofs - Creating initial metabolite set


This repository is for all Crossomics Metabolite set creation scripts. The main scripts will be run from run_array_job.sh (creating the initial sets) and run.sh (extending the sets).



INITIAL METABOLITE SET CREATION
GLOBAL: loop: run_array_job.sh 
		--> array_job.sh 
			--> findSurroundings.R 
				<--> getMetsPathwayCommons.R
					<--> getLeftOrRight.R
			save Gene.RData in folder 'mss'
	Repeat

* run_array_job.sh loops through numbers corresponding to the amount of genes present (see point 4). It calls on array_job.sh and gives array_job.sh a SGE_TASK_ID.

* array_job.sh runs the R script 'findSurroundings.R' and gives it three values: 1. the SGE_TASK_ID as gene_id (a single number), 2. an output folder and 3. a folder where all supporting scripts can be found. 

* R Scripts that are used by / downstream of findSurroundings.R are getMetsPathwayCommons.R and getLeftOrRight.R

* The number of genes that is looped through is dependent on how many have been downloaded from the HGNC biomart website (https://biomart.genenames.org/, accessed 5 june 2019). All genes have been downloaded with the attributes 'HGNC data: Approved symbol' and 'Protein resources: Uniprot accession', leading to 43200 genes. These have been filtered so both attributes have a value (many genes miss a Uniprot accession value), leading to a total of 20351 genes left. These 20351 genes is the number that is looped over by run_array_job.sh. The list of genes is also saved and loaded in by findSurroundings.R

* For every gene, a connection with the Pathway Commons website is made to extract gene-metabolite information. If no metabolites can be found for a gene, the gene is skipped and no output file is made.

* The output of this process is a folder (mss) which contains the initial metabolite sets. 


POSSIBLE IMPROVEMENTS FOR FUTURE
This process has been run locally on my MAC. It would have been best to perform this pipeline on the hpc, however that is not possible at the moment as rJava can't be loaded inside the qsub/qlogin environment.

Save gene.RData files in gene.RDS files. The RDS-format is a tad easier to work with, although performance should stay the same.




EXTENDING METABOLITE SETS
run.sh calls on extend.sh and gives it a SGE_TASK_ID, which in turn calls extend.R and gives it a gene_id (a single number, the same as the SGE_TASK_ID), an output folder and a folder where all supporting scripts can be found.
Input used by extend.R are the files created by run_array_job.sh (mss_WG_step_0) and Recon2.RData.
Scripts used by / downstream of extend.R are findMetabolicEnvironmentLocal.R, getMets2ReconID.R, convert.R, getReactionsRecon.R, removeMetsFromSet.R, getPreviousNext.R, and getMetsRecon2V2.R.
Output: 4 different folders with extended metabolite sets (0 to 4 steps away from the primary reaction).
It is unknown whether this can be run on the hpc, but probably not as it probably depends on rJava somewhere.