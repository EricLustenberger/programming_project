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

% Create grid of unemployment benefits for the analysis
mu_min = 0.01;
mu_max = 0.6;
mu_n = 20;
mu = linspace(mu_min, mu_max, mu_n);
mu(7) = mu(6); % the original point does not converge
mu(6) = 0.15; % the original point does not converge
mu(8) = 0.18; % the original point does not converge

% Corresponding PI_UE grid for analysis 
PI_UE_grid = [0.418006431, 0.413867753, 0.409823146, 0.405844156, 0.401954115,0.4, 0.398125746, 0.396341463, 0.38708909, 0.383537395, 0.380061395, 0.376636922, 0.373284328, 0.369980363, 0.366744717, 0.363555009, 0.360430298, 0.35734902, 0.354329636, 0.351351351];