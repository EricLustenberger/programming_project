%% Create graph of welfare measures for different values of unemployment benefit (mu) with mu = 0.15 as baseline
% Produces the welfare vs benefits plots recalling the results obtained
% from the welfare_analysis.m file. 

% Add path to matlab-json parser
addpath ../library/matlab-json/
json.startup

% Load unemployment insurance grid
baseline = json.read(project_paths('IN_MODEL_SPECS', ['baseline.json']));

mu_grid = baseline.mu_grid;
mu_n = length(mu_grid);

% Load the data
for i=1:mu_n
    filename = ['baseline_' num2str(i) '.mat'];
    c(i) = load(project_paths('OUT_ANALYSIS',filename), 'c');
    k(i) = load(project_paths('OUT_ANALYSIS',filename), 'k');
    wk(i) = load(project_paths('OUT_ANALYSIS',filename), 'wk');
    wc(i) = load(project_paths('OUT_ANALYSIS',filename), 'wc');
    c_equiv(i,:) = [c(i).c.equivalent_mean, c(i).c.equivalent_median, c(i).c.equivalent_unemployed_mean, c(i).c.equivalent_unemployed_median, c(i).c.equivalent_employed_mean, c(i).c.equivalent_employed_median];
    k_equiv(i,:) = [k(i).k.equivalent_mean, k(i).k.equivalent_median, k(i).k.equivalent_unemployed_mean, k(i).k.equivalent_unemployed_median, k(i).k.equivalent_employed_mean, k(i).k.equivalent_employed_median];
end

output_baseline = 3.3539;
% Obtain cash equivalent relative to baseline output
rel_k_equiv = k_equiv./output_baseline; % Get cash equivalent relative to output
      
% Create the figures
figure (1)
plot(mu_grid(1:17)', k_equiv(1:17,3)./ output_baseline,'m', mu_grid(1:17)', k_equiv(1:17,5)./ output_baseline,'g' ...
    ,mu_grid(1:17)', k_equiv(1:17,4)./ output_baseline,'--m', mu_grid(1:17)', k_equiv(1:17,6)./ output_baseline,'--g')
legend('unemployed','employed')
xlabel('unemployment benefit')
ylabel('cash equivalent / output')
refline (0,0)
axis tight
saveas(gcf,project_paths('OUT_FIGURES','cash_equivalent_ue_vs_e.png'));

figure (2)
plot(mu_grid(1:13)', k_equiv(1:13,1)./ output_baseline,'r', mu_grid(1:13)', k_equiv(1:13,2)./ output_baseline,'--r')
legend('mean','median')
xlabel('unemployment benefit')
ylabel('cash equivalent / output')
refline (0,0)
saveas(gcf,project_paths('OUT_FIGURES','cash_equivalent_total.png'));

figure (3)
plot(mu_grid', c_equiv(:,3),'m', mu_grid', c_equiv(:,5),'g' ...
    ,mu_grid', c_equiv(:,4),'--m', mu_grid', c_equiv(:,6),'--g')
legend('unemployed','employed')
xlabel('unemployment benefit')
ylabel('consumption equivalent')
refline (0,1)
axis tight
saveas(gcf,project_paths('OUT_FIGURES','consumption_equivalent_ue_vs_e.png'));

figure (4)
plot(mu_grid', c_equiv(:,1),'r', mu_grid', c_equiv(:,2),'--r')
legend('mean','median')
xlabel('unemployment benefit')
ylabel('consumption equivalent')
refline (0,1)
axis tight
saveas(gcf,project_paths('OUT_FIGURES','consumption_equivalent_total.png'));

figure (5)
plot(mu_grid(1:9)', k_equiv(1:9,1)./ output_baseline,'r', mu_grid(1:9)', (c_equiv(1:9,1)-1).*100, 'g' ...
    ,mu_grid(1:9)', k_equiv(1:9,2)./ output_baseline,'--r' , mu_grid(1:9)', (c_equiv(1:9,2)-1).*100,'--g')
legend('cash equivalent', 'consumption equivivalent')
xlabel('unemployment benefit')
%ylabel('cash equivalent / output', 'consumption equivalent')
refline (0,0)
axis tight
saveas(gcf,project_paths('OUT_FIGURES','consumption_vs_cash_equivalent.png'));


for i= baseline.wealth_plotting
% Wealth vs equivalent figures 
% consumption equivalent
	figure (5+i)
	plot(wc(i).wc.k_emp(wc(i).wc.sort_index_emp),c(i).c.equivalent_emp_sorted,'g.',wc(i).wc.k_unemp(wc(i).wc.sort_index_unemp),c(i).c.equivalent_unemp_sorted,'r.')
	legend('employed','unemployed')
	xlabel('wealth')
	ylabel('consumption equivalent')
	refline (0,1)
	saveas(gcf,project_paths('OUT_FIGURES','consumption_equivalent_vs_wealth.png'));

	figure (6+i)
	plot(wk(i).wk.k_emp(wk(i).wk.sort_index_emp),k(i).k.equivalent_emp_sorted./output_baseline...
    	,'g.',wk(i).wk.k_unemp(wk(i).wk.sort_index_unemp),k(i).k.equivalent_unemp_sorted./output_baseline,'r.')
	legend('employed','unemployed')
	xlabel('wealth')
	ylabel('cash equivalent / output')
	refline (0,0)
	saveas(gcf,project_paths('OUT_FIGURES','cash_equivalent_vs_wealth.png'));
end




