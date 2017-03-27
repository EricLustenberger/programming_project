%% Welfare analysis of two Aiyagari models with different parameters
% Main file for the welfare analysis. Calculates means and medians of the
% consumption- and cash equivalent and saves them to be recalled by the
% graph_equivalent.m file for the walfare vs benefit plots. The means and
% medians are both calculated separately for the unemployed and employed as
% well as in the total. 
% Moreover, when not looping over the different values for mu and PI_UE,
% this file prints the welfare vs wealth plots for a predetermined policy
% change. 

% of the unemployed and employed, as well as the total means and medians welfare analysis,

clear
close all


% Add path to relevant model code
addpath ../model_code/

% Add path to matlab-json parser
addpath ../library/matlab-json/
json.startup

% Load model specifications
baseline = json.read(project_paths('IN_MODEL_SPECS', ['baseline.json']));
par.beta = baseline.discount_factor;
par.sigma = baseline.risk_aversion; 
par.mu = baseline.baseline_policy;
par.k_min= baseline.borrowing_constraint;
par.delta = baseline.depreciation_rate;
par.alpha = baseline.output_elasticity_of_capital;
par.z = baseline.productivity; 
par.PI_EU = baseline.job_loss_probability;
par.PI_UE = baseline.job_finding_probability;
mu = baseline.replacement_rate;

% solution methods
% Default setup for analysis
method.HH = 'FP'; % 'FP' for fixed-point iteration in household problem, 'FPend' for using endogenous grid points
method.sim = 'simulation'; % 'simulation' for simulation, 'histogram' for histogram method
method.agg = 'bisectio'; % 'bisection' to use bisection method, gradual updating otherwise

setup % load setup

% Corresponding PI_UE grid for analysis 
PI_UE_grid = [0.418006431, 0.413867753, 0.409823146, 0.405844156, 0.401954115,0.4, 0.398125746, 0.396341463, 0.38708909, 0.383537395, 0.380061395, 0.376636922, 0.373284328, 0.369980363, 0.366744717, 0.363555009, 0.360430298, 0.35734902, 0.354329636, 0.351351351];

% If you want to compare one steady state with the baseline, change
% parameters for i=2 and leave i=1 as the baseline.

% Get steady state utility and welfare measures for the model with two different parameters
% i=1, calculates the expected livetime utility of the baseline
% i=2, calculates the expeced livetime utility of a one-time permanent change in the benefit level for
% the values specified in the mu-grid
for i=1:2
    if i==1 
        [ k.one, c.one, K.one, sim.one, store.one, mat.one, grid.one ] = aiyagari_solver( par, func, method );
        U.one.guess = func.U(c.one.guess);
        U.one.lifetime = zeros(2,grid.one.k_no); %expected life time utility
        dist=100;
        while dist>1e-8
            EU = par.PI * U.one.lifetime; 
            Unew = U.one.guess' + par.beta * EU;
            dist = max(max(abs(Unew - U.one.lifetime)));
            U.one.lifetime = Unew;
        end
        U.one.lifetime(U.one.lifetime == -Inf) = -999999; % Get rid of -Inf for negative consumption for extrapolation
        ind_no = size(sim.one.k,2);
        T = size(sim.one.k,1);
        U.one.extrap = NaN(ceil(T/2),ind_no);
        for t = ceil((T+1)/2):T % Extrapolate life time utility
            U.one.extrap(t-ceil(T/2),sim.one.e(t,:)==1) = interp1(grid.one.k, U.one.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
            U.one.extrap(t-ceil(T/2),sim.one.e(t,:)==2) = interp1(grid.one.k, U.one.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');
        end
    elseif i==2 
        for ii=1:size(mu,2) 
            tic
            iteration = ii
            par.mu = mu(ii);
            %par.PI_UE = PI_UE_grid(ii); %loop over the job-finding
            %probability grid for the experiment

        % adjust method for different values of mu to ensure convergence
            if ii < 9
                method.HH = 'FP'; % Depending on the mu, you might have to change it to 'FP' or 'FPend' to converge;
                method.agg = 'bisectio'; % Depending on the mu, you might have to change it to 'bisection' or 'bisectio' to converge;
            else
                method.HH = 'FP'; % Depending on the mu, you might have to change it to 'FP' or 'FPend' to converge;
                method.agg = 'bisection'; % Depending on the mu, you might have to change it to 'bisection' or 'bisectio' to converge;
            end
            setup % refresh setup for new parameter
            [ k.two, c.two, K.two, sim.two, store.two, mat.two, grid.two ] = aiyagari_solver( par, func, method );
            U.two.guess = func.U(c.two.guess);
            U.two.lifetime = zeros(2,grid.two.k_no); %expected life time utility
            dist=100;
            while dist>1e-8
                EU = par.PI * U.two.lifetime; 
                Unew = U.two.guess' + par.beta * EU;
                dist = max(max(abs(Unew - U.two.lifetime)));
                U.two.lifetime = Unew;
            end
            U.two.lifetime(U.two.lifetime == -Inf) = -999999; % Get rid of -Inf for extrapolation
            ind_no = size(sim.one.k,2);
            T = size(sim.one.k,1);
            U.two.extrap = NaN(ceil(T/2),ind_no);
            for t = ceil((T+1)/2):T % Extrapolate life time utility with idiosyncratic transition
                U.two.extrap(t-ceil(T/2),sim.one.e(t,:)==1) = interp1(grid.one.k, U.two.lifetime(1,:), sim.one.k(t,sim.one.e(t,:)==1), 'linear', 'extrap');
                U.two.extrap(t-ceil(T/2),sim.one.e(t,:)==2) = interp1(grid.one.k, U.two.lifetime(2,:), sim.one.k(t,sim.one.e(t,:)==2), 'linear', 'extrap');
            end
            % Calculate the consumption equivalentof the policy reform. If c.equivalent > 1,
            % agents prefer the policy change. If 1 > c.equivalent > 0, agents prefer the
            % baseline model.
            [ c, wc ] = consumption_equivalent (par, method, sim, U); %calling consumption_equivalent fct.

            % Calculate the cash equivalent of the policy change, If
            % k.equivalent > 0, agents prefer the policy change. If
            % k.equivalent < 0, agents prefer the baseline model.
            [ k, wk ] = cash_equivalent (method, grid, sim, U); %calling cash_equivalent fct.

            % Construct a variable of values to be stored
            keep.c.equivalent_mean = c.equivalent_mean;
            keep.c.equivalent_median = c.equivalent_median;
            keep.c.equivalent_unemployed_mean = c.equivalent_unemployed_mean;
            keep.c.equivalent_unemployed_median = c.equivalent_unemployed_median;
            keep.c.equivalent_employed_mean = c.equivalent_employed_mean;
            keep.c.equivalent_employed_median = c.equivalent_employed_median;
            keep.c.equivalent_emp_sorted = c.equivalent_emp_sorted;
            keep.c.equivalent_unemp_sorted = c.equivalent_unemp_sorted;
            keep.k.equivalent_mean = k.equivalent_mean;
            keep.k.equivalent_median = k.equivalent_median;
            keep.k.equivalent_unemployed_mean = k.equivalent_unemployed_mean;
            keep.k.equivalent_unemployed_median = k.equivalent_unemployed_median;
            keep.k.equivalent_employed_mean = k.equivalent_employed_mean;
            keep.k.equivalent_employed_median = k.equivalent_employed_median;
            keep.k.equivalent_emp_sorted = k.equivalent_emp_sorted;
            keep.k.equivalent_unemp_sorted = k.equivalent_unemp_sorted;
            keep.K = K.two.guess;
            keep.wc.k_unemp = wc.k_unemp;
            keep.wc.k_emp = wc.k_emp;
            keep.wc.sort_index_unemp = wc.sort_index_unemp;
            keep.wc.sort_index_emp = wc.sort_index_emp;
            keep.wk.k_unemp = wk.k_unemp;
            keep.wk.k_emp = wk.k_emp;
            keep.wk.sort_index_unemp = wk.sort_index_unemp;
            keep.wk.sort_index_emp = wk.sort_index_emp;


            filename = ['baseline_' num2str(ii) '.mat']; %Store the values for each grid-point
            save(project_paths('OUT_ANALYSIS', filename),'-struct','keep');  
           
            toc
        end
    end
end
