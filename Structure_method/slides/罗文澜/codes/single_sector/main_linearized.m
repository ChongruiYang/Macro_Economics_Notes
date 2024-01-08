% This scritp solves the main model using log linearized methods
clear;

% Add path
addpath('functions');

% Set parameters
params = setup;

% Set data
data = get_data;

% Modify the migration matrix of RoW to be strongly connected
N = data.N;
real_small = 1e-5;
mu_ni_m10 = data.mu_ni_m1;

mu_ni_m1 = mu_ni_m10;
mu_ni_m1_domestic = mu_ni_m1(1:N-1,1:N-1);
mu_ni_m1(1:N-1,N) = real_small;
mu_ni_m1_domestic = mu_ni_m1_domestic - real_small/(N-1);
mu_ni_m1(1:N-1,1:N-1) = mu_ni_m1_domestic;
mu_ni_m1(N,1:N-1) = real_small / (N-1);
mu_ni_m1(N,N) = 1 - real_small;
data.mu_ni_m1 = mu_ni_m1;

% Solve the linearized model
eq = solve_linearized(params,data);

% Apply shocks to solve benchmark and counter-factuals
% Baseline
idx = eq.idx;
T = data.T;
numShocks = size(eq.RR,2);
ln_A_n_t = cumsum(log(data.dot_A_n_t),2);
ln_kappa_ni_t = cumsum(log(data.dot_kappa_ni_t),3);
ln_tau_ni_t = cumsum(log(data.dot_tau_ni_t),3);
ln_kappa_in_n_t = reshape(sum(eq.SS.*ln_kappa_ni_t, 2),[N,T]);
ln_kappa_out_n_t = reshape(sum(eq.TT.*permute(ln_kappa_ni_t,[2,1,3]), 2),[N,T]);
ln_tau_in_n_t = reshape(sum(eq.EE.*permute(ln_tau_ni_t,[2,1,3]), 2),[N,T]);
ln_tau_out_n_t = reshape(sum(eq.DD.*ln_tau_ni_t, 2),[N,T]);

ln_shocks_t = zeros(numShocks,T);
ln_shocks_t(idx.A,:) = ln_A_n_t;
ln_shocks_t(idx.kappa_in,:) = ln_kappa_in_n_t;
ln_shocks_t(idx.kappa_out,:) = ln_kappa_out_n_t;
ln_shocks_t(idx.tau_in,:) = ln_tau_in_n_t;
ln_shocks_t(idx.tau_out,:) = ln_tau_out_n_t;
ln_L_n_t = zeros(N,T);
for t=2:T-1
    ln_L_n_t(:,t+1) = eq.HH*ln_L_n_t(:,t) + eq.RR*ln_shocks_t(:,t);
end
baseline.ln_L_n_t = ln_L_n_t;

% Counterfactual
ln_A_n_t = cumsum(log(data.hat_A_n_t.*data.dot_A_n_t),2);
ln_kappa_ni_t = cumsum(log(data.hat_kappa_ni_t.*data.dot_kappa_ni_t),3);
ln_tau_ni_t = cumsum(log(data.hat_tau_ni_t.*data.dot_tau_ni_t),3);
ln_kappa_in_n_t = reshape(sum(eq.SS.*ln_kappa_ni_t, 2),[N,T]);
ln_kappa_out_n_t = reshape(sum(eq.TT.*permute(ln_kappa_ni_t,[2,1,3]), 2),[N,T]);
ln_tau_in_n_t = reshape(sum(eq.EE.*permute(ln_tau_ni_t,[2,1,3]), 2),[N,T]);
ln_tau_out_n_t = reshape(sum(eq.DD.*ln_tau_ni_t, 2),[N,T]);

ln_shocks_t = zeros(numShocks,T);
ln_shocks_t(idx.A,:) = ln_A_n_t;
ln_shocks_t(idx.kappa_in,:) = ln_kappa_in_n_t;
ln_shocks_t(idx.kappa_out,:) = ln_kappa_out_n_t;
ln_shocks_t(idx.tau_in,:) = ln_tau_in_n_t;
ln_shocks_t(idx.tau_out,:) = ln_tau_out_n_t;
ln_L_n_t = zeros(N,T);
for t=2:T-1
    ln_L_n_t(:,t+1) = eq.HH*ln_L_n_t(:,t) + eq.RR*ln_shocks_t(:,t);
end
cf.ln_L_n_t = ln_L_n_t;

save('mat/linearized.mat','baseline','cf','data','eq');