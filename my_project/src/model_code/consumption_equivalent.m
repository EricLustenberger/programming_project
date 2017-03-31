function [ c, wc ] = consumption_equivalent (par, sim, U)
%{
CONSUMPTION EQUIVALENT is a function that calculates the consumption equivalent of a policy change.
%}  

    % Calculate the consumption equivalent
    % If value > 1, agents prefer steady state one, if 1 > value > 0, agents prefer
    % policy change ( = steady state two).
    % Consumption equivalent tested against model with policy change:
    c.equivalent = ((U.two.extrap.*(1-par.sigma).*(1-par.beta)+1)...
        ./(U.one.extrap.*(1-par.sigma).*(1-par.beta)+1)).^(1/(1-par.sigma));

    % Calculate the mean and median consumption equivalent
    c.equivalent_mean = mean(mean(c.equivalent));
    c.equivalent_median = median(median(c.equivalent));
        
    % Get consumption equivalent for employed and unemployed
    T = size(sim.one.e,1);
    for t=1:ceil(T/2)
        c.equivalent_unemployed_mean(t,:) = mean(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
        c.equivalent_unemployed_median(t,:) = median(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
        c.equivalent_employed_mean(t,:) = mean(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
        c.equivalent_employed_median(t,:) = median(c.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
    end
    c.equivalent_unemployed_mean = mean(c.equivalent_unemployed_mean);
    c.equivalent_unemployed_median = median(c.equivalent_unemployed_median);
    c.equivalent_employed_mean = mean(c.equivalent_employed_mean);
    c.equivalent_employed_median = median(c.equivalent_employed_median);
        
    % Sort the consumption equivalent for employed and unemployed in the last period to plot them
    % against wealth
    [c.equivalent_unemp_sorted, wc.sort_index_unemp] = sort(c.equivalent(end,sim.one.e(T,:)==1),'descend');
    [c.equivalent_emp_sorted, wc.sort_index_emp] = sort(c.equivalent(end,sim.one.e(T,:)==2),'descend');
    wc.k_unemp = sim.one.k(T,sim.one.e(T,:)==1);
    wc.k_emp = sim.one.k(T,sim.one.e(T,:)==2);
               
end

