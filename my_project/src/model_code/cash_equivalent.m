function [ k, wk ] = cash_equivalent( grid, sim, U )
%{
CASH EQUIVALENT is a function that calculates the capital equivalent of a policy change.
%}   
    ind_no = size(sim.one.k,2);
    T = size(sim.one.k,1);
    k.equivalent = NaN(ceil(T/2),ind_no);
    k.compensated = NaN(ceil(T/2),ind_no);
        
    % Extrapolate the cash equivalent
    for t = 1:ceil(T/2)
        k.compensated(t,sim.one.e(t+ceil(T/2),:)==1) = interp1(U.one.lifetime(1,:), grid.one.k, U.two.extrap(t,sim.one.e(t+ceil(T/2),:)==1), 'linear', 'extrap');
        k.compensated(t,sim.one.e(t+ceil(T/2),:)==2) = interp1(U.one.lifetime(2,:), grid.one.k, U.two.extrap(t,sim.one.e(t+ceil(T/2),:)==2), 'linear', 'extrap');
        k.equivalent(t,:) = k.compensated(t,:) - sim.one.k(t+ceil(T/2),:);   
        k.equivalent_unemployed_mean(t,:) = mean(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
        k.equivalent_unemployed_median(t,:) = median(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==1));
        k.equivalent_employed_mean(t,:) = mean(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
        k.equivalent_employed_median(t,:) = median(k.equivalent(t,sim.one.e(ceil(T/2)+t,:)==2));
    end
        
    % Calculate mean and median cash equivalent
    k.equivalent_mean = mean(mean(k.equivalent));
    k.equivalent_median = median(median(k.equivalent));
        
    % Get cash equivalent for employed and unemployed
    k.equivalent_unemployed_mean = mean(k.equivalent_unemployed_mean);
    k.equivalent_unemployed_median = median(k.equivalent_unemployed_median);
    k.equivalent_employed_mean = mean(k.equivalent_employed_mean);
    k.equivalent_employed_median = median(k.equivalent_employed_median);
        
    % Sort the cash equivalent for employed and unemployed in the last period to plot them against
    % wealth
    [k.equivalent_unemp_sorted, wk.sort_index_unemp] = sort(k.equivalent(end,sim.one.e(T,:)==1),'descend');
    [k.equivalent_emp_sorted, wk.sort_index_emp] = sort(k.equivalent(end,sim.one.e(T,:)==2),'descend');
    wk.k_unemp = sim.one.k(T,sim.one.e(T,:)==1);
    wk.k_emp = sim.one.k(T,sim.one.e(T,:)==2);
           
end

