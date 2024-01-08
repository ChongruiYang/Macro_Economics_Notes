% This script solves the baseline economy of a one-sector model with
% roundabout production and dynamic migration decision
% Author: Wenlan Luo
% Date: 7/7/2021
clear;

% Add path
addpath('functions');

% Set parameters
params = setup;

% Set data
data = get_data;

% 
params.TOL_DYNAMIC_EQ = 1e-5;
eq = solve_trans(params,data);

save('mat/baseline.mat','eq','data','params');
