%% Set parameters
par.beta = 0.99; % discount factor
par.sigma = 2; % risk aversion
par.mu = 0.15; % replacement rate of unemployed
par.k_min = 1e-15; % borrowing constraint as share of capital stock

% firms (production function F(K,L)=z*K^alpha*L^(1-alpha)
par.delta = 0.025; % depreciation rate
par.alpha = 0.36; % output elasticity of capital
par.z = 1; % productivity

% transition probabilities
par.PI_UE = 0.4; % baseline job-finding probability
par.PI_EU = 0.044444444; % fixed layoff probability at baseline level

% solution methods
% Default setup for analysis
method.HH = 'FP'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
method.sim = 'simulation'; % 'simulation' for simulation, 'histogram' for histogram method
method.agg = 'bisectio'; % 'bisection' to use bisection method, gradual updating otherwise