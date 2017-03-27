function [ k, c, K, sim, store, mat, grid] = aiyagari_solver( par, func, method, simulation)
% AIYAGARI MODEL: Heterogeneous agents model due to idiosyncratic labour
% shocks. Agents self-sinsure against unemploment by building capital
% stock.
% Input variables:
%   par = (calibrated) Parameters that decribe the economy.
%   func = Helpful functions to automate certain calculations.
%   method = Describe the method of iteration/ simutation.

    %redefine simulation
    sim.k = simulation.k;
    sim.e = simulation.e;
    sim.shock = simulation.shock;
    ind_no = simulation.ind_no;
    T = simulation.T;

    K.rep = func.K(1/par.beta-1+par.delta);% capital stock of representative agent useful for comparison

    grid.k_no = 100; % number of grid points for agents' capital holdings
    grid.k = linspace(par.k_min*K.rep,3*K.rep,grid.k_no); % grid for agents' policy function

    grid.dist_no = 1000;
    grid.dist = linspace(grid.k(1),grid.k(end),grid.dist_no); % grid for distribution of agents - use finer grid

    % useful matrices
    mat.k = [grid.k',grid.k']; % replicate grid for unemployed and employed
    mat.income = @(K) func.w(K)*repmat([par.mu,1-par.tau],grid.k_no,1); % matrix with income of each agent

    K.lims = [K.rep,grid.k(end)]; % initial limits for bisection method

    % initial guesses
    k.guess = mat.k; % policy function of agents (1st column unemployed, 2nd column employed)
    if strcmp(method.agg,'bisection')
        K.guess = (K.lims(1)+K.lims(2))/2;
    else
        K.guess = K.rep; % aggregate capital stock
    end

    d1 = 1;
    iter = 0;
    while d1>1e-6 && iter<50 % loop for aggregate problem
        iter = iter+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%    Solve household problem given prices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d2 = 1;
        while d2>1e-10 % loop for household problem

            if strcmp(method.HH,'FP') % Fixed-point iteration     
                % calculate capital chosen next period, depending on current and
                % future employment status
                for i=1:2 % employment this period
                    for j=1:2 % employment next period
                        % capital choice next period from policy function
                        k.next(i,:,j) = interp1(grid.k,k.guess(:,j),k.guess(:,i),'linear','extrap'); 
                    end
                        % consumption next period from budget constraint
                        c.next(i,:,:) = max(1e-10,(1+func.r(K.guess)-par.delta)*[k.guess(:,i),k.guess(:,i)] - squeeze(k.next(i,:,:)) + mat.income(K.guess));
                end

                % calculate expected marginal utility of consumption next period
                Emuc_next(:,1) = par.PI(1,1)*func.muc(c.next(1,:,1))+par.PI(1,2)*func.muc(c.next(1,:,2));  % currently unemployed
                Emuc_next(:,2) = par.PI(2,1)*func.muc(c.next(2,:,1))+par.PI(2,2)*func.muc(c.next(2,:,2));  % currently employed       

                % calculate implied consumption this period from Euler equation
                c.current = func.muc_inv(par.beta*(1+func.r(K.guess)-par.delta)*Emuc_next);

                % calculate implied capital demand from budget constraint
                k.new = (1+func.r(K.guess)-par.delta)*mat.k + mat.income(K.guess) - c.current;

            elseif strcmp(method.HH,'FPend') % Fixed-point iteration with endogenous grid points method
                % calculate consumption next period, if agents are on the grid in
                % the next period
                c.next = max(1e-10,(1+func.r(K.guess)-par.delta)*mat.k-k.guess+mat.income(K.guess));

                % calculate expected marginal utility of consumption next period,
                % if future capital holdings are on the grid
                Emuc_next(:,1) = par.PI(1,1)*func.muc(c.next(:,1))+par.PI(1,2)*func.muc(c.next(:,2));  % currently unemployed
                Emuc_next(:,2) = par.PI(2,1)*func.muc(c.next(:,1))+par.PI(2,2)*func.muc(c.next(:,2));  % currently employed

                % calculate implied consumption this period from Euler equation
                c.current = func.muc_inv(par.beta*(1+func.r(K.guess)-par.delta)*Emuc_next);

                % calculate implied capital stock this period from budget constraint
                k.current = (c.current+mat.k-mat.income(K.guess))./(1+func.r(K.guess)-par.delta);

                % invert policy function to get k_next(k_current) on the grid for k
                for i=1:2 % currently unemployed and employed
                    k.new(:,i) = interp1(k.current(:,i),grid.k,grid.k,'linear','extrap');
                end
            end

            % apply borrowing constraint to get new policy function
            k.new = max(par.k_min*K.guess,k.new); 

            d2 = norm(abs(k.new-k.guess)./(1+abs(k.guess))); % deviation between guess and new policy function

            % update policy function
            k.guess = k.guess + 0.5*(k.new-k.guess);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%    Find distribution of agents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
         
        sim.e(1,1:round(par.L*ind_no))=2; % initial individuals that are employed
        sim.k(1,:) = K.guess; % initial capital holdings

        for t=2:T
            sim.e(t,sim.e(t-1,:)==1) = 1+(sim.shock(t,sim.e(t-1,:)==1)<=par.PI(1,2)); % new employment status of previously unemployed
            sim.e(t,sim.e(t-1,:)==2) = 1+(sim.shock(t,sim.e(t-1,:)==2)<=par.PI(2,2)); % new employment status of previously employed    
            sim.k(t,sim.e(t,:)==1) = interp1(grid.k,k.guess(:,1),sim.k(t-1,sim.e(t,:)==1),'linear','extrap'); % capital demand of currently unemployed
            sim.k(t,sim.e(t,:)==2) = interp1(grid.k,k.guess(:,2),sim.k(t-1,sim.e(t,:)==2),'linear','extrap'); % capital demand of currently employed
        end

            K.demand = mean(mean(sim.k(ceil(T/2):end,:))); % average capital holdings over second half of sample
            sim.L = mean(mean(sim.e(ceil(T/2):end,:)==2)); % average employment
        

        d1 = abs(K.demand-K.guess)./(1+K.guess); % deviation between guess for capital and demanded capital stock

        store.K_guess(iter)=K.guess;
        store.K_next(iter)=K.demand;

        if K.demand>K.guess % update limits of bisection interval
            K.lims(1) = K.guess;
        else
            K.lims(2) = K.guess;
        end

        % update guess for aggregate capital stock
        if strcmp(method.agg,'bisection')
            K.guess = (K.lims(1)+K.lims(2))/2;
        else
            K.guess = K.guess + 0.05*(K.demand-K.guess);
        end
        disp(['Iteration: ',num2str(iter),', K_guess: ',num2str(K.guess),', K_demand: ',num2str(K.demand)])

        c.guess = (1+func.r(K.guess)-par.delta)*mat.k-k.guess+mat.income(K.guess); % Get consumption policy function
        c.guess = max(0, c.guess);
    end


    disp('_____________________________________________________________')
    disp('Aggregate variables (log-deviation from representative agent)')
    disp(['Capital:     ',num2str(K.guess),' (',num2str(log(K.guess/K.rep)),')'])
    disp(['Output:      ',num2str(func.Y(K.guess)),' (',num2str(log(func.Y(K.guess)/func.Y(K.rep))),')'])
    disp(['Consumption: ',num2str(func.C(K.guess)),' (',num2str(log(func.C(K.guess)/func.C(K.rep))),')'])

end

