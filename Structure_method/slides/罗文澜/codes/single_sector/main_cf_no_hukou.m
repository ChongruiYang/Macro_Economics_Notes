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
data = get_data;
% Shutdown the change in trade costs
data.hat_kappa_ni_t(:) = 1;

% Load baseline
baseline = load('mat/baseline.mat');

% 
eq = solve_cf_trans(baseline.eq,params,data);

save('mat/cf_no_hukou.mat','eq','params','data');