% This script exports additional statistics of the counterfactual
% equilibrium
clear;
baseline = load('mat/baseline.mat');
cf = load('mat/cf.mat');

% Construct population dynamics of the baseline at year 10 and year 20
pop_change_baseline_10 = baseline.eq.L_n_t(:,6) ./ baseline.eq.L_n_t(:,1) - 1;
pop_change_baseline_20 = baseline.eq.L_n_t(:,21) ./ baseline.eq.L_n_t(:,1) - 1;

% Construct population change at the impact period, year 5, year 10 and
% year 20
outStruct.pop_change_cf_to_baseline_10 = (1 - cf.eq.prime_L_n_t(1:end-1,11) ./ baseline.eq.L_n_t(1:end-1,11))*100;
outStruct.pop_change_cf_to_baseline_20 = (1 - cf.eq.prime_L_n_t(1:end-1,21) ./ baseline.eq.L_n_t(1:end-1,21))*100;
outStruct.citygb = baseline.data.cityList;

writetable(struct2table(outStruct),'tables/cf_pop_change_by_city.csv')