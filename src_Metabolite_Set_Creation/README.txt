README - Crossomics Metabolite Set Creation

31-05-2019 - Marten Kerkhofs



This repository is for all Crossomics Metabolite set creation scripts. The main scripts will be run from run_array_job.sh (creating the initial sets) and run.sh (extending the sets).


run_array_job.sh calls on array_job.sh and gives it a SGE_TASK_ID, which in turn calls findSurroundings.R and gives it a gene_id (a single number, the same as the SGE_TASK_ID), an output folder and a folder where all supporting scripts can be found. 
Scripts that are used by / downstream of findSurroundings.R are getMetsPathwayCommons.R and getLeftOrRight.R
The output of this process is a folder (mss_WG_step_0) which contains the initial metabolite sets. 
Although this function should be run on the hpc, that is not possible at the moment as rJava can't be loaded inside the qsub/qlogin environment.


run.sh calls on extend.sh and gives it a SGE_TASK_ID, which in turn calls extend.R and gives it a gene_id (a single number, the same as the SGE_TASK_ID), an output folder and a folder where all supporting scripts can be found.
Input used by extend.R are the files created by run_array_job.sh (mss_WG_step_0) and Recon2.RData.
Scripts used by / downstream of extend.R are findMetabolicEnvironmentLocal.R, getMets2ReconID.R, convert.R, getReactionsRecon.R, removeMetsFromSet.R, getPreviousNext.R, and getMetsRecon2V2.R.
Output: 4 different folders with extended metabolite sets (0 to 4 steps away from the primary reaction).
It is unknown whether this can be run on the hpc, but probably not as it probably depends on rJava somewhere.