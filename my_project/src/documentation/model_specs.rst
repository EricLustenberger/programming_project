.. _model_specifications:


Model specifications
======================

The model specifications display all exogenously given parameters. The directory **src.model_specs** contains `JSON <http://www.json.org/>`_ files with model specifications. The choice of JSON is motivated by the attempt to be language-agnostic: JSON is quite expressive and there are parsers for nearly all languages. The two files, baseline.json and endogenous_job_finding_duration.json differ with respect to the job finding probability. The contains a grid of the job finding probability, thus the when looping over the different policy changes, the job finding probability will also adopt - unlike in the baseline case, where PI_UE is set at its baseline value. 

Explaining the Parameters
-------------------------
 
    * *beta* - the discount factor 
    * *sigma* - the risk aversion parameter
    * *k_min* - the minimum capital each household must hold in each period
    * *mu* - the baseline replacement rate 
    * *delta* - the discount factor
    * *alpha* - the output elasticity of capital
    * *z* - the productivity
    * *PI_EU* - fixed layoff probability at baseline 
    * *PI_UE* - the baseline job-finding probability or unemployment duration
    * *mu_grid* - are the different unemployment benefit policies considered 
    * *wealth_plotting* - allows to determine which policy change should be plotted 
    * *k_no* - the grid for capital holdings
    * *dist_no* - the grid for the distribution of agents
    * *ind_no* - number of agents considered 
    * *T* - number of time periods considered
    * *PI_UE_grid* - for the analysis, the grid to let the unemployment duration fluctuate 

