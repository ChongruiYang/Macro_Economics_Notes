function eq = solve_cf_trans(baseline,params,data)
% SOLVE_CF_TRANS solve the transition path of the counterfactual economy

% Unpack from parameters and data
DAMP_P = params.DAMP_P;
DAMP_U = params.DAMP_U;
TOL_DYNAMIC_EQ = params.TOL_DYNAMIC_EQ;
TOL_TEMPORARY_EQ = params.TOL_TEMPORARY_EQ;
PRINT_FREQ = params.PRINT_FREQ;

beta = params.beta;
nu = params.nu;
theta = params.theta;

gamma_n = data.gamma_n;
gamma_tilde_n = data.gamma_tilde_n;

N = data.N;
T = data.T;

hat_A_n_t = data.hat_A_n_t;
hat_kappa_ni_t = data.hat_kappa_ni_t;
hat_tau_ni_t = data.hat_tau_ni_t;

% Inititate the path
prime_pi_ni_t = ones(N,N,T);
prime_X_n_t = ones(N,T);
prime_mu_ni_t = ones(N,N,T);
prime_L_n_t = ones(N,T);

prime_pi_ni_t(:,:,1) = baseline.pi_ni_t(:,:,1);
prime_X_n_t(:,1) = baseline.X_n_t(:,1);
prime_mu_ni_t(:,:,1) = baseline.mu_ni_t(:,:,1);
prime_L_n_t(:,1) = baseline.L_n_t(:,1);

hat_u_n_t = ones(N,T);
hat_L_n_t = ones(N,T);
hat_P_n_t = ones(N,T);
hat_w_n_t = ones(N,T);

iter = 0;
converged = false;
while ~converged
    iter = iter+1;
    % Solve the change in migration share
    for t=1:T
        % Migration share unchanged at impact to reflect that the shocks are
        % unexpected
        if t==1
            hat_u_n_t(:,1) = 1;
            prime_mu_ni_t(:,:,1) = baseline.mu_ni_t(:,:,1);
            continue;
        end

        % Extract from time series
        prime_mu_ni_last = prime_mu_ni_t(:,:,t-1);
        dot_mu_ni = baseline.dot_mu_ni_t(:,:,t);
        dot_mu_ni(isnan(dot_mu_ni)) = 1;    % RoW
        hat_tau_ni = hat_tau_ni_t(:,:,t);
        
        if t<T
            hat_u_n_next = hat_u_n_t(:,t+1);
        else
            hat_u_n_next = ones(N,1); % assume T+1 at the steady state
        end
        
        if t==2
            mu1_tilde = baseline.mu_ni_t(:,:,2).*reshape(hat_u_n_t(:,2).^(beta/nu),[1,N]);
            numer = mu1_tilde.*reshape(hat_u_n_next.^(beta/nu),[1,N]) .* hat_tau_ni.^(-1/nu);
            denom = sum(numer,2);
            prime_mu_ni_t(:,:,2) = numer ./ denom;
            continue;
        end
       
        numer = prime_mu_ni_last .* dot_mu_ni .* reshape(hat_u_n_next.^(beta/nu),[1,N]) .* hat_tau_ni.^(-1/nu);
        denom = sum(numer, 2);
        prime_mu_ni = numer ./ denom;

        % Assign to time series variables
        prime_mu_ni_t(:,:,t) = prime_mu_ni;
    end
    % clear temporary variables to avoid mistakes
    clear prime_mu_ni_last dot_mu_ni hat_tau_ni prime_mu_ni
    
    % Solve the change in labor
    for t=1:T
        if t<=2
            prime_L_n_t(:,t) = baseline.L_n_t(:,t);
            hat_L_n_t(:,t) = 1;
        else
            % Extract from time series
            prime_mu_ni_last = prime_mu_ni_t(:,:,t-1);
            prime_L_n_last = prime_L_n_t(:,t-1);
            
            % Move forward
            prime_L_n = permute( sum(prime_mu_ni_last.*reshape(prime_L_n_last,[N,1]), 1), [2,1] );
            
            % Assign to time series variables
            prime_L_n_t(:,t) = prime_L_n;
            hat_L_n_t(:,t) = (prime_L_n./prime_L_n_last) ./ baseline.dot_L_n_t(:,t);
        end
    end
    clear prime_mu_ni_last prime_L_n_last prime_L_n;
    
    % Solve the temporary equilibrium forward
    for t=2:T
        % Extract from time series
        dot_pi_ni = baseline.dot_pi_ni_t(:,:,t);
        dot_w_n = baseline.dot_w_n_t(:,t);
        dot_L_n = baseline.dot_L_n_t(:,t);
        
        hat_kappa_ni = hat_kappa_ni_t(:,:,t);
        hat_A_n = hat_A_n_t(:,t);
        
        hat_L_n = hat_L_n_t(:,t);
        prime_pi_ni_last = prime_pi_ni_t(:,:,t-1);
        prime_X_n_last = prime_X_n_t(:,t-1);
        
        % initialize
        hat_w_n = hat_w_n_t(:,t);
        hat_P_n = hat_P_n_t(:,t);

        % Share from the last period
        prime_wL_n_last = gamma_n .* permute(sum(prime_pi_ni_last.*prime_X_n_last,1),[2,1]);
        
        % solve the temoprary equilibrium in changes
        minor_converged = false;
        while ~minor_converged
            % Construct the change
            hat_x_n = hat_w_n.^gamma_n .* hat_P_n.^gamma_tilde_n;
            numer = prime_pi_ni_last.*dot_pi_ni .* (reshape(hat_x_n,[1,N]).*hat_kappa_ni).^(-theta) ...
                .* reshape(hat_A_n.^(theta*gamma_n),[1,N]);
            sum_numer = sum(numer, 2);
            hat_P_n_new = sum_numer.^(-1/theta);
            prime_pi_ni = numer ./ sum_numer;
            
            % Form the linear coefs to solve X
            prime_pi_in = prime_pi_ni;
            linearCoefs_ni = eye(N) - ( reshape(gamma_n,[N,1]).* permute(prime_pi_in,[2,1]) ...
                + reshape(gamma_tilde_n,[N,1]).*permute(prime_pi_in,[2,1]) );
            % Normalization
            linearCoefs_ni(N+1,:) = 1;
            cons_n = zeros(N+1,1);
            cons_n(N+1) = 1;
            % Solve
            prime_X_n = linearCoefs_ni \ cons_n;
            
            % Calculate implied wage
            hat_w_n_new = gamma_n.* permute(sum( prime_pi_ni.*prime_X_n, 1 ), [2,1]) ...
                ./ prime_wL_n_last ./ (dot_w_n.*dot_L_n) ./ hat_L_n;
            
            % Calculate metrics
            metric_w = max(abs(log(hat_w_n_new(:)) - log(hat_w_n(:))));
            metric_P = max(abs(log(hat_P_n_new(:)) - log(hat_P_n(:))));
            metric = max([metric_w,metric_P]);
            
            minor_converged = metric<TOL_TEMPORARY_EQ;
            
            % Update
            hat_w_n = exp(DAMP_P*log(hat_w_n_new) + (1-DAMP_P)*log(hat_w_n));
            hat_P_n = exp(DAMP_P*log(hat_P_n_new) + (1-DAMP_P)*log(hat_P_n));
        end
        
        prime_pi_ni_t(:,:,t) = prime_pi_ni;
        prime_X_n_t(:,t) = prime_X_n;
        hat_w_n_t(:,t) = hat_w_n;
        hat_P_n_t(:,t) = hat_P_n;
        fprintf('.');
    end
    fprintf('\n');
    
    % Some calculations before updating
    hat_omega_n_t = hat_w_n_t ./ hat_P_n_t;
    
    % Update dot_u backward
    hat_u_n_t_new = ones(N,T);
    for t=T:-1:1
        % Utility unchanged before impact to reflect that the shocks are
        % unexpected
        if t==1
            hat_u_n_t_new(:,1) = 1;
            continue;
        end

        % Extract from time series
        prime_mu_ni_last = prime_mu_ni_t(:,:,t-1);
        dot_mu_ni = baseline.dot_mu_ni_t(:,:,t);
        dot_mu_ni(isnan(dot_mu_ni)) = 1;    % RoW
        hat_tau_ni = hat_tau_ni_t(:,:,t);
        hat_omega_n = hat_omega_n_t(:,t);
        
        if t<T
            hat_u_n_next = hat_u_n_t(:,t+1);
        else
            hat_u_n_next = ones(N,1); % assume T+1 at the steady state
        end
        
        % Allow for jump of u at the impact
        if t==2
            mu1_tilde = baseline.mu_ni_t(:,:,2).*reshape(hat_u_n_t(:,2).^(beta/nu),[1,N]);
            numer = mu1_tilde .* reshape(hat_u_n_next.^(beta/nu),[1,N]) .* hat_tau_ni.^(-1/nu);
        else
            % Update backward
            numer = prime_mu_ni_last.*dot_mu_ni .* reshape(hat_u_n_next.^(beta/nu),[1,N]) .* hat_tau_ni.^(-1/nu);
        end
        hat_u_n_new = hat_omega_n .* sum( numer, 2 ).^nu;
        
        % Assign to time series
        hat_u_n_t_new(:,t) = hat_u_n_new;
    end
    clear prime_mu_ni_last dot_mu_ni hat_tau_ni hat_omega_n hat_u_n_new
    
    % Metric
    metric_dot_u = max(abs(log(hat_u_n_t(:)) - log(hat_u_n_t_new(:))));
    metric = metric_dot_u;
    
    converged = metric < TOL_DYNAMIC_EQ;
    
    % Update
    hat_u_n_t = DAMP_U*hat_u_n_t_new + (1-DAMP_U)*hat_u_n_t;
    
    % Print
    if mod(iter,PRINT_FREQ)==0 || converged
        fprintf('Major iter, metric: %d, %g\n', iter, metric);
    end
end

eq = v2struct(prime_pi_ni_t,prime_X_n_t,prime_L_n_t,prime_mu_ni_t,...
    hat_u_n_t,hat_L_n_t,hat_P_n_t,hat_w_n_t,hat_omega_n_t);

end