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

% Load baseline
baseline = load('mat/baseline.mat');

% 
eq = solve_cf_trans(baseline.eq,params,data);

save('mat/cf.mat','eq','params','data');