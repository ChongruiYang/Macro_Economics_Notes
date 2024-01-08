function eq = solve_trans(params,data)
% SOLVE_TRANS solve the transition path of the baseline economy
% Author: Wenlan Luo; luowl@sem.tsinghua.edu.cn
% Date: 7/7/2021

% Unpack from parameters and data
DAMP_U = params.DAMP_U;
TOL_DYNAMIC_EQ = params.TOL_DYNAMIC_EQ;
PRINT_FREQ = params.PRINT_FREQ;

beta = params.beta;
nu = params.nu;

N = data.N;
T = data.T;

dot_A_n_t = data.dot_A_n_t;
dot_kappa_ni_t = data.dot_kappa_ni_t;
dot_tau_ni_t = data.dot_tau_ni_t;

% Inititate the path
pi_ni_t = ones(N,N,T);
X_n_t = ones(N,T);
L_n_t = ones(N,T);
mu_ni_t = ones(N,N,T);

pi_ni_t(:,:,1) = data.pi_ni_0;
X_n_t(:,1) = data.X_n_0;
L_n_t(:,1) = data.L_n_0;
mu_ni_m1 = data.mu_ni_m1;

dot_u_n_t = ones(N,T);
dot_L_n_t = ones(N,T);
dot_P_n_t = ones(N,T);
dot_w_n_t = ones(N,T);
dot_mu_ni_t = ones(N,N,T);
dot_pi_ni_t = ones(N,N,T);

iter = 0;
converged = false;
while ~converged
    iter = iter+1;
    % Solve the migration share implied by dot_u forward
    for t=1:T
        % Extract from time series
        if t==1
            mu_ni_last = mu_ni_m1;
        else
            mu_ni_last = mu_ni_t(:,:,t-1);
        end
        
        if t<T
            dot_u_n_next = dot_u_n_t(:,t+1);
        else
            dot_u_n_next = ones(N,1); % assume T+1 at the steady state
        end
        
        dot_tau_ni = dot_tau_ni_t(:,:,t);
        
        mu_ni_numerator = mu_ni_last .* reshape(dot_u_n_next.^(beta/nu),[1,N]) .* dot_tau_ni.^(-1/nu);
        mu_ni_denominator = sum(mu_ni_numerator, 2);
        mu_ni = mu_ni_numerator ./ mu_ni_denominator;

        % Assign to time series variables
        mu_ni_t(:,:,t) = mu_ni;
        dot_mu_ni_t(:,:,t) = mu_ni ./ mu_ni_last;
    end
    % clear temporary variables to avoid mistakes
    clear mu_ni_last dot_u_n_next dot_tau_ni mu_ni;
    
    % Solve the implied labor forward
    for t=2:T
        % Extract from time series
        mu_ni_last = mu_ni_t(:,:,t-1);
        L_n_last = L_n_t(:,t-1);
        
        % Move forward
        L_n = permute( sum(mu_ni_last.*reshape(L_n_last,[N,1]), 1), [2,1] );
        dot_L_n = L_n ./ L_n_last;
        
        % Assign to time series variables
        L_n_t(:,t) = L_n;
        dot_L_n_t(:,t) = dot_L_n;
    end
    clear mu_ni_last L_n_last L_n dot_L_n;
    
    % Solve the temporary equilibrium forward
    for t=2:T
        % Extract from time series
        pi_ni_last = pi_ni_t(:,:,t-1);
        X_n_last = X_n_t(:,t-1);
        dot_L_n = dot_L_n_t(:,t);
        dot_A_n = dot_A_n_t(:,t);
        dot_kappa_ni = dot_kappa_ni_t(:,:,t);
        
        % Use solutions from the last iteration as initial guesses
        initial.dot_P_n = dot_P_n_t(:,t);
        initial.dot_w_n = dot_w_n_t(:,t);
        
        % Solve forward for one time step
        [dot_P_n,dot_w_n,pi_ni,X_n] = solve_temporary_eq(pi_ni_last,X_n_last,dot_L_n,dot_A_n,dot_kappa_ni,params,data,initial);
        
        % Assign to time series
        dot_P_n_t(:,t) = dot_P_n;
        dot_w_n_t(:,t) = dot_w_n;
        pi_ni_t(:,:,t) = pi_ni;
        X_n_t(:,t) = X_n;
        dot_pi_ni_t(:,:,t) = pi_ni ./ pi_ni_last;
        fprintf('.');
    end
    fprintf('\n');
    clear pi_ni_last X_n_last dot_L_n dot_A_n dot_kappa_ni;
    
    % Some calculations before updating
    dot_omega_n_t = dot_w_n_t ./ dot_P_n_t;
    
    % Update dot_u backward
    dot_u_n_t_new = ones(N,T);
    for t=T:-1:1
        % Extract from time series
        if t==1
            mu_ni_last = mu_ni_m1;
        else
            mu_ni_last = mu_ni_t(:,:,t-1);
        end
        
        if t<T
            dot_u_n_next = dot_u_n_t(:,t+1);
        else
            dot_u_n_next = ones(N,1); % assume T+1 at the steady state
        end
        
        dot_omega_n = dot_omega_n_t(:,t);
        dot_tau_ni = dot_tau_ni_t(:,:,t);
        
        % Update backward
        dot_u_n_new = dot_omega_n .* sum( mu_ni_last .* reshape(dot_u_n_next.^(beta/nu),[1,N]) .* dot_tau_ni.^(-1/nu), 2 ).^nu;
        
        % Assign to time series
        dot_u_n_t_new(:,t) = dot_u_n_new;
    end
    clear mu_ni_last dot_u_n_next dot_omega_n dot_tau_ni dot_u_n_new
    
    % Metric
    metric_dot_u = max(abs(log(dot_u_n_t(:)) - log(dot_u_n_t_new(:))));
    metric = metric_dot_u;
    
    converged = metric < TOL_DYNAMIC_EQ;
    
    % Update
    dot_u_n_t = DAMP_U*dot_u_n_t_new + (1-DAMP_U)*dot_u_n_t;
    
    % Print
    if mod(iter,PRINT_FREQ)==0 || converged
        fprintf('Major iter, metric: %d, %g\n', iter, metric);
    end
end

eq = v2struct(pi_ni_t,X_n_t,L_n_t,mu_ni_t,...
    dot_u_n_t,dot_L_n_t,dot_P_n_t,dot_w_n_t,dot_omega_n_t,dot_mu_ni_t,dot_pi_ni_t);

end