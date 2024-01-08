% This script solves the counter factual of the baseline economy applying
% the dynamic hat algebra
% Author: Wenlan Luo
% Date: 7/7/2021

% Add path
close all;
addpath('functions');

% Set parameters
params = setup;

% Set data
data = get_data_test;

% Load baseline
baseline = solve_trans(params,data);
cf = solve_cf_trans(baseline,params,data);

%%%%%%%%%% Labor dynamics 
figure;
subplot(2,2,2);
plot(cf.prime_L_n_t','LineWidth',1.5);
title('Labor, counterfactual');
xlabel('period');
ylabel('level');
ylimCurrent = ylim;
legend({'Region 1','Region 2'});

subplot(2,2,1);
plot(baseline.L_n_t');
title('Labor, baseline');
xlabel('period');
ylabel('level');
ylim(ylimCurrent);
% print('figures/test_cf_pop.png','-dpng');

%%%%%% Example to convert dot variable to level relative to period 0
baseline.omega_n_t = cumprod(baseline.dot_omega_n_t,2);
cf.omega_n_t = cumprod(baseline.dot_omega_n_t.*cf.hat_omega_n_t,2);
% figure;
subplot(2,2,4);
plot(cf.omega_n_t','LineWidth',1.5);
title('Real wage, counterfactual');
xlabel('period');
ylabel('level to period 0');
ylimCurrent = ylim;
legend({'Region 1','Region 2'});

subplot(2,2,3);
plot(baseline.omega_n_t');
title('Real wage, baseline');
xlabel('period');
ylabel('level to period 0');
ylim(ylimCurrent);
print('figures/test_cf_pop_real_wage.png','-dpng');
