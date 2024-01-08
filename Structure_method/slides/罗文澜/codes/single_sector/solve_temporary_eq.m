function [dot_P_n,dot_w_n,pi_ni,X_n] = solve_temporary_eq(pi_ni_last,X_n_last,dot_L_n,dot_A_n,dot_kappa_ni,params,data,initial)
% SOLVE_TEMPORARY_EQ solves the temporary equilibrium given informaion
% required
% Author: Wenlan Luo; luowl@sem.tsinghua.edu.cn
% Date: 7/7/2021

% Extract params
DAMP_P = params.DAMP_P;
TOL_TEMPORARY_EQ = params.TOL_TEMPORARY_EQ;

theta = params.theta;

N = data.N;

gamma_n = data.gamma_n;
gamma_tilde_n = data.gamma_tilde_n;

% Use Warm up externally provided
dot_P_n = initial.dot_P_n;
dot_w_n = initial.dot_w_n;

% Evaluate equilibrium objects from the last step
wL_n_last = gamma_n .* permute(sum(pi_ni_last.*X_n_last,1),[2,1]);

converged = false;
iter = 0;
while ~converged
    iter = iter + 1;
    
    % Construct the implied trade share
    dot_x_n = dot_w_n.^gamma_n .* dot_P_n.^gamma_tilde_n;
    pi_ni_numerator = pi_ni_last .* (reshape(dot_x_n,[1,N]).*dot_kappa_ni).^(-theta) ...
        .* reshape(dot_A_n.^(theta*gamma_n),[1,N]);
    sum_pi_ni_numerator = sum(pi_ni_numerator, 2);
    dot_P_n_new = sum_pi_ni_numerator.^(-1/theta);
    pi_ni = pi_ni_numerator ./ sum_pi_ni_numerator;
    
    % Form the linear coefs to solve X
    pi_in = pi_ni;
    linearCoefs_ni = eye(N) - ( reshape(gamma_n,[N,1]).* permute(pi_in,[2,1]) ...
        + reshape(gamma_tilde_n,[N,1]).*permute(pi_in,[2,1]) );
    % Normalization
    linearCoefs_ni(N+1,:) = 1;    
    cons_n = zeros(N+1,1);
    cons_n(N+1) = 1;
    % Solve
    X_n = linearCoefs_ni \ cons_n;
    
    % Calculate implied wage
    dot_w_n_new = gamma_n.* permute(sum( pi_ni.*X_n, 1 ), [2,1]) ...
        ./ wL_n_last ./ dot_L_n;
    
    % Calculate metrics
    metric_w = max(abs(log(dot_w_n_new(:)) - log(dot_w_n(:))));
    metric_P = max(abs(log(dot_P_n_new(:)) - log(dot_P_n(:))));
    metric = max([metric_w,metric_P]);
    
    converged = metric<TOL_TEMPORARY_EQ;
    
    % Update
    dot_w_n = exp(DAMP_P*log(dot_w_n_new) + (1-DAMP_P)*log(dot_w_n));
    dot_P_n = exp(DAMP_P*log(dot_P_n_new) + (1-DAMP_P)*log(dot_P_n));
end
end