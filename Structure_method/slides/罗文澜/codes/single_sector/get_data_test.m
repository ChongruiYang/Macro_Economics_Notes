function data = get_data(params)
% GET_DATA returns data for a test problem
T = 100;
N = 10;

%%%%%%%%%% Input-output tables
gamma_n = 0.5*ones(N,1);
gamma_tilde_n = 0.5*ones(N,1);

%%%%%%%%%% Initial shares
% import share 20%
damp = 0.2 / (N-1);
pi_ni_0 = damp*ones(N,N);
pi_ni_0(eye(N,'logical')) = 1 - damp*(N-1);

% annual migration rate 5%
damp = 0.05 / (N-1);
mu_ni_m1 = damp*ones(N,N);
mu_ni_m1(eye(N,'logical')) = 1 - damp*(N-1);

% symmetric region
X_n_0 = ones(N,1) / N;
L_n_0 = ones(N,1) / N;

%%%%%%%%%%%%%% Shocks
% baseline
dot_A_n_t = ones(N,T);
dot_kappa_ni_t = ones(N,N,T);
dot_tau_ni_t = ones(N,N,T);

% counterfactuals
hat_A_n_t = ones(N,T);
hat_kappa_ni_t = ones(N,N,T);
hat_tau_ni_t = ones(N,N,T);
hat_A_n_t(1,2) = 0.9;

clear params;
data = v2struct;
end