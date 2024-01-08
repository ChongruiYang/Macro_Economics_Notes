function eq = solve_linearized(params,data)
% SOLVE_LINEARIZED solve the log linearized coefficients

% Unpack from parameters and data
theta = params.theta;
beta = params.beta;
nu = params.nu;

N = data.N;
T = data.T;
gamma_n = data.gamma_n;
gamma_tilde_n = data.gamma_tilde_n;

% Only one time shock is considered
ln_A_n_t = cumsum(log(data.hat_A_n_t),2);
% ln_A_n_t = ln_A_n_t - ln_A_n_t(:,end);
ln_kappa_ni_t = cumsum(log(data.hat_kappa_ni_t),3);
ln_tau_ni_t = cumsum(log(data.hat_tau_ni_t),3);

% Inititate the path
pi_ni_0 = data.pi_ni_0;
X_n_0 = data.X_n_0;
L_n_0 = data.L_n_0;
mu_ni_m1 = data.mu_ni_m1;


%%%%%%%%%%%%%%%%%%%%%%%%% Construct Share Matrix %%%%%%%%%%%%%%%%%%%%

% Construct the share matrix, based on trade flows and population flows at
% the initial period
% All matrix has row sum of one (i.e., sum(xx,2)==1)
% Do tranposition carefully

% Expenditure share
SS = pi_ni_0;
X_ni_0 = X_n_0 .* pi_ni_0;  % X_ni defined as trade flows i->n
% Income share
TT = permute(X_ni_0 ./ sum(X_ni_0, 1), [2,1]);
% Out migration share
DD = mu_ni_m1;
L_ni_0 = L_n_0 .* mu_ni_m1; % L_ni defined as population flows n->i
% In migration share
EE = permute(L_ni_0 ./ sum(L_ni_0, 1), [2,1]);

%
X_n_share_0 = X_n_0 / sum(X_n_0);
L_n_share_0 = L_n_0 / sum(L_n_0);

% Average trade and migration costs
ln_kappa_in_n_t = reshape(sum(SS.*ln_kappa_ni_t, 2),[N,T]);
ln_kappa_out_n_t = reshape(sum(TT.*permute(ln_kappa_ni_t,[2,1,3]), 2),[N,T]);
ln_tau_in_n_t = reshape(sum(EE.*permute(ln_tau_ni_t,[2,1,3]), 2),[N,T]);
ln_tau_out_n_t = reshape(sum(DD.*ln_tau_ni_t, 2),[N,T]);

%%%%%%%%%%%%%%%%%%%%%%%% Static problem %%%%%%%%%%%%%%%%%%%%%%%%%

% Shock index in the column vector f
idx = struct;

ptr = 0;

idx.A = ptr+1:ptr+N;
ptr = ptr+N;

idx.kappa_in = ptr+1:ptr+N;
ptr = ptr+N;

idx.kappa_out = ptr+1:ptr+N;
ptr = ptr+N;

idx.tau_in = ptr+1:ptr+N;
ptr = ptr+N;

idx.tau_out = ptr+1:ptr+N;
ptr = ptr+N;

numShocks = ptr;


% Variable index for temporary equilibrium
ptr = 0;

idx.x = ptr+1:ptr+N;
ptr = ptr+N;

idx.P = ptr+1:ptr+N;
ptr = ptr+N;

idx.X = ptr+1:ptr+N;
ptr = ptr+N;

idx.w = ptr+1:ptr+N;
ptr = ptr+N;

numTemporaryVars = ptr;

% Form the linear coefficients for temporary equilibrium
Omega = zeros(numTemporaryVars,numTemporaryVars);
Lambda = zeros(numTemporaryVars,N);
Xi = zeros(numTemporaryVars,numShocks);

% Equation by equation
Omega(idx.x, idx.x) = eye(N);
Omega(idx.x, idx.w) = -diag(gamma_n);
Omega(idx.x, idx.P) = -diag(gamma_tilde_n);

Omega(idx.P, idx.P) = eye(N);
Omega(idx.P, idx.x) = -SS;
Xi(idx.P, idx.A) = -SS*diag(gamma_n);
Xi(idx.P, idx.kappa_in) = eye(N);

Omega(idx.X, idx.X) = eye(N) - TT;
Omega(idx.X, idx.x) = theta*eye(N);
Omega(idx.X, idx.P) = -theta*TT;
Xi(idx.X, idx.A) = theta*diag(gamma_n);
Xi(idx.X, idx.kappa_out) = -theta*eye(N);

Omega(idx.w, idx.w) = eye(N);
Omega(idx.w, idx.X) = -eye(N);
Lambda(idx.w, :) = -eye(N);

% One normalize equation. We normalize sum(X_n) = 1;
idx_normalize = size(Omega,1) + 1;
Omega(idx_normalize,:) = 0;
Omega(idx_normalize,idx.X) = X_n_share_0(:)';
Xi(idx_normalize,:) = 0;
Lambda(idx_normalize,:) = 0;
if (rank(Omega)~=4*N)
    error('Rank error for temporary equilibrium');
end

% Solve for coefficients for the contemporary equilibrium
BB = Omega \ Lambda;
CC = Omega \ Xi;
BB_real_wage = BB(idx.w,:) - BB(idx.P,:);
CC_real_wage = CC(idx.w,:) - CC(idx.P,:);

%%%%%%%%%%%%%%%%%%%%%% Dynamic problem %%%%%%%%%%%%%%%%%%%%%%%%%
ptr = 0;

idx.L = ptr+1:ptr+N;
ptr = ptr+N;

idx.u = ptr+1:ptr+N;
ptr = ptr+N;

Gamma = zeros(2*N,2*N);
Gamma(idx.L,idx.L) = eye(N);
Gamma(idx.L,idx.u) = -beta/nu*(eye(N)-EE*DD);
Gamma(idx.u,idx.u) = beta*DD;

Theta = zeros(2*N,2*N);
Theta(idx.L,idx.L) = EE;
Theta(idx.u,idx.u) = eye(N);
Theta(idx.u,idx.L) = -BB_real_wage;

Psi = zeros(2*N,numShocks);
Psi(idx.L,idx.tau_in) = -1/nu*eye(N);
Psi(idx.L,idx.tau_out) = 1/nu*EE;
Psi(idx.u,:) = -CC_real_wage;
Psi(idx.u,idx.tau_out) = Psi(idx.u,idx.tau_out) + eye(N);

% Normalization
Gamma_ex = Gamma;
normalize_idx = 2*N+1;
Gamma_ex(normalize_idx,idx.L) = L_n_share_0';
Theta_ex = Theta;
Theta_ex(normalize_idx,idx.L) = 0;
Psi_ex = Psi;
Psi_ex(normalize_idx,:) = 0;
Theta_tilde = Gamma_ex \ Theta_ex;
Psi_tilde = Gamma_ex \ Psi_ex;

% %{
% Check the consistent of solution given by solab
[FF,HH] = solab(eye(2*N), Theta_tilde, N);
RR_u = (eye(N) + FF*Theta_tilde(idx.L,idx.u) - Theta_tilde(idx.u,idx.u)) \ (Psi_tilde(idx.u,:) - FF*Psi_tilde(idx.L,:));
RR = Theta_tilde(idx.L,idx.u)*RR_u + Psi_tilde(idx.L,:);

check_consis1 = Theta_tilde(idx.L,idx.u)*FF + Theta_tilde(idx.L,idx.L) - HH;
check_consis2 = Theta_tilde(idx.u,idx.L) + Theta_tilde(idx.u,idx.u)*FF - FF*HH;
%}

% %{
% Use gensys
% [HH,CC,RR,fmat,fwt,ywt,gev,eu,loose] = gensys(Gamma,Theta,Psi(:,1),zeros(2*N,1),zeros(2*N,1));
% [HH,CC,RR,fmat,fwt,ywt,gev,eu,loose] = gensys(eye(2*N),Theta_tilde,Psi_tilde(:,1),zeros(2*N,1),zeros(2*N,1));
% [HH,~,RR,fmat,fwt,ywt,gev,eu,loose] = gensys(eye(2*N),Theta_tilde,zeros(2*N,1),Psi_tilde,zeros(2*N,2*N));
% [HH,~,RR,fmat,fwt,ywt,gev,eu,loose] = gensys(eye(2*N),Theta_tilde,zeros(2*N,1),zeros(2*N,1),zeros(2*N,1));
% [HH_full,~,RR,fmat,fwt,ywt,gev,eu,loose] = gensys(Gamma,Theta,zeros(2*N,1),zeros(2*N,1),zeros(2*N,1));
% [HH,~] = gensysToAMA(eye(2*N),Theta_tilde,zeros(2*N,1),zeros(2*N,1),zeros(2*N,1),[],'both');
% [HH,~,RR,fmat,fwt,ywt,gev,eu,loose] = gensys(Gamma_ex,Theta_ex,zeros(2*N,1),Psi_ex,zeros(2*N,2*N));
% RR_future = ywt*inv(I-fmat*inv(L))*fwt;
%}

%{
HH = HH_full(1:N,1:N);
FF = 
RR_u = (eye(N) + FF*Theta_tilde(idx.L,idx.u) - Theta_tilde(idx.u,idx.u)) \ (Psi_tilde(idx.u,:) - FF*Psi_tilde(idx.L,:));
RR = Theta_tilde(idx.L,idx.u)*RR_u + Psi_tilde(idx.L,:);
%}


%{
% We stack the shock also into states and use solab to solve the difference
% equation
idx.L_ex = numShocks + idx.L;
idx.u_ex = numShocks + idx.u;
idx.shock_ex = 1:numShocks;
Gamma_expanded = zeros(2*N+numShocks,2*N+numShocks);
Gamma_expanded([idx.L_ex,idx.u_ex],[idx.L_ex,idx.u_ex]) = eye(2*N);
Gamma_expanded(idx.shock_ex,idx.shock_ex) = eye(numShocks);

Theta_expanded = zeros(2*N+numShocks,2*N+numShocks);
Theta_expanded(idx.shock_ex,idx.shock_ex) = eye(numShocks);
Theta_expanded([idx.L_ex,idx.u_ex],1:numShocks) = Psi_tilde;
Theta_expanded([idx.L_ex,idx.u_ex],[idx.L_ex,idx.u_ex]) = Theta_tilde;
[FF,PP] = solab(Gamma_expanded, Theta_expanded, numShocks+N);

% Extract the state transition function
HH = PP(idx.L_ex,idx.L_ex);
RR = PP(idx.L_ex,idx.shock_ex);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate forward %%%%%%%%%%%%%%%%%
% Stacking shocks
ln_shocks_t = zeros(numShocks,T);
ln_shocks_t(idx.A,:) = ln_A_n_t;
ln_shocks_t(idx.kappa_in,:) = ln_kappa_in_n_t;
ln_shocks_t(idx.kappa_out,:) = ln_kappa_out_n_t;
ln_shocks_t(idx.tau_in,:) = ln_tau_in_n_t;
ln_shocks_t(idx.tau_out,:) = ln_tau_out_n_t;

% %{
ln_L_n_t = zeros(N,T);
% ln_L_n_t(:,1) = data.ln_L_n_0;
for t=2:T-1
    ln_L_n_t(:,t+1) = HH*ln_L_n_t(:,t) + RR*ln_shocks_t(:,t);
end
%}

%{
% Un comment the following for results from gensys
ln_vars_t = zeros(2*N,T);
for t=2:T
%     ln_vars_t(:,t) = HH*ln_vars_t(:,t-1) + RR*ln_shocks_t(:,t);
    ln_vars_t(:,t) = HH*ln_vars_t(:,t-1) + CC*ln_A_n_t(1,t);
end
ln_L_n_t = ln_vars_t(idx.L,:);
%}

eq = v2struct(SS,TT,EE,DD,BB,CC,HH,RR,ln_L_n_t,BB_real_wage,CC_real_wage,idx,ln_shocks_t);

end