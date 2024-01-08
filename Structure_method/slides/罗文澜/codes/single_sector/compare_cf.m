% Compare results
clear;
% main_trial;
test_baseline;
test_cf;
test_linearized;
baseline = load('mat_test/baseline.mat');
nonlinear = load('mat_test/cf.mat');
linearized = load('mat_test/linearized.mat');

% linearized.ln_omega_n_t = linearized.eq.BB_real_wage*linearized.ln_L_n_t + linearized.CC_real_wage*linearized.ln_shocks_t;

figure; hold on;
plot((log(nonlinear.eq.prime_L_n_t) - log(nonlinear.eq.prime_L_n_t(:,1)))','k');
plot((linearized.eq.ln_L_n_t+log(baseline.eq.L_n_t) - log(baseline.eq.L_n_t(:,1)))','r');
