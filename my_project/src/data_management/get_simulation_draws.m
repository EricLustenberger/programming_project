%{ 
Draw simulated samples from two uncorrelated uniform variables 
(locations in two dimensions) for two types of agents and store 
them in a 3-dimensional matrix. Note: In principle, one would read 
the number of dimensions etc. from the "IN_MODEL_SPECS" file, this is 
to demonstrate the most basic use of run_m_script only. 
%}

simulation.ind_no = 20; % number of individuals simulated
simulation.T = 20; % number of periods simulated

rng('default') % reset random number generator

simulation.k = zeros(simulation.T,simulation.ind_no); % simulated values of capital stock
simulation.e = ones(simulation.T,simulation.ind_no); % simulated employment status
simulation.shock = rand(simulation.T,simulation.ind_no); % shocks for employment transition

save(project_paths('OUT_DATA','simulation.mat'),'simulation')