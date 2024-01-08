% This script solves the linearized model and compares with the nonlinear
% solutions using dynamic hat algebra
% Author: Wenlan Luo; luowl@sem.tsinghua.edu.cn
% Date: 7/7/2021

clear;
close all;
addpath('functions');

% Set parameters
params = setup;
params.PRINT_FREQ = 20;

% Set data
data = get_data_test;

% Solve different versions
params.TOL_DYNAMIC_EQ = 1e-5;
data.hat_A_n_t(1,2) = 0.5;
baseline = solve_trans(params,data);
cf = solve_cf_trans(baseline,params,data);
linearized = solve_linearized(params,data);

% Compare
figure; hold on;
plot((log(cf.prime_L_n_t(1,:)) - log(baseline.L_n_t(1,1)))');
plot(linearized.ln_L_n_t(1,:)','--');
legend({'Nonlinear Solution','Log Linearized Solution'});
xlabel('period');
title('Change in log(population) of region 1');
print('figures/test_linearized_shock_large.png','-dpng');

%%%%%%%%%%%%%%%%% Set a smaller shock %%%%%%%%%
params.TOL_DYNAMIC_EQ = 1e-5;
rng(1123);
data.hat_A_n_t(:,2) = ones(10,1);
data.hat_A_n_t(1,2) = 0.9;
% data.hat_A_n_t(:,2) = 1-rand(10,1)*0.05;
% data.hat_tau_ni_t(:,:,2) = rand(10,10);
% data.hat_tau_ni_t(:,:,2) = ones(10,10);
% data.hat_kappa_ni_t(:,:,2) = 1-rand(10,10)*0.1;
baseline = solve_trans(params,data);
cf = solve_cf_trans(baseline,params,data);
linearized = solve_linearized(params,data);

figure; hold on;
plot((log(cf.prime_L_n_t(1,:)) - log(baseline.L_n_t(1,:)))');
plot(linearized.ln_L_n_t(1,:)','--');
legend({'Nonlinear Solution','Log Linearized Solution'});
xlabel('period');
title('Change in log(population) of region 1');
print('figures/test_linearized_shock_small.png','-dpng');




