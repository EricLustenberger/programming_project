%{ 
Draw simulated samples from two uncorrelated uniform variables 
(locations in two dimensions) for two types of agents and store 
them in a 3-dimensional matrix. Note: In principle, one would read 
the number of dimensions etc. from the "IN_MODEL_SPECS" file, this is 
to demonstrate the most basic use of run_m_script only. 
%}

ind_no = 20; % number of individuals simulated
T = 20; % number of periods simulated

rng('default') % reset random number generator
sim.k = zeros(T,ind_no); % simulated values of capital stock
sim.e = ones(T,ind_no); % simulated employment status
sim.shock = rand(T,ind_no); % shocks for employment transition

n_types = 2;
n_draws = 30000;

% Set random seed
rng(12345, 'twister')

% 3d matrix with random draws from a (0,1) uniform distribution
sample = rand(n_draws, 2, n_types);

% Save 
save(project_paths('OUT_DATA', 'sample.mat'), 'sample');
