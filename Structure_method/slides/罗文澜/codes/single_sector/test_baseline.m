% This script solves the baseline economy of a one-sector model with
% roundabout production and dynamic migration decision
% Author: Wenlan Luo; luowl@sem.tsinghua.edu.cn
% Date: 7/7/2021
clear;

% Add path
addpath('functions');

% Set parameters
params = setup;

% Set data
data = get_data_test;
eq = solve_trans(params,data);

% Inspect the dynamics of labor
figure;
plot(eq.L_n_t');
title('Labor dynamics');
xlabel('period');
ylabel('level');

% Set some anticipated shock, so initial period is not at steady state
data.dot_A_n_t(1,2) = 0.9;  % Lower region 1's productivity by 10%
eq_new = solve_trans(params,data);
figure;
plot(eq_new.L_n_t','LineWidth',1.5);
title('Labor');
xlabel('period');
ylabel('level');
legend({'Region 1','Region 2'});
print('figures/test_baseline_pop.png','-dpng');

%%%%%% Example to convert dot variable to level relative to period 0
omega_n_t = cumprod(eq_new.dot_omega_n_t,2);
figure;
plot(omega_n_t','LineWidth',1.5);
title('Real wage');
xlabel('period');
ylabel('level relative to period 0');
legend({'Region 1','Region 2'},'Location','Best','FontSize',12);
print('figures/test_baseline_real_wage.png','-dpng');