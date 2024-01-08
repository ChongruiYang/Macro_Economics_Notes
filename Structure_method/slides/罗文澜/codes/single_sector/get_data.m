function data = get_data_JT(params)
%%%%%%%%%%%%% Get data from Fan (2019) and Fan, Lu and Luo (2021)
% Trade flows
tradeMat = load('data_Luo/bilateral_trade_flows.mat');

% Migration flows
migration_cityid = load('data_JT/cityid.mat');
migration_no_reform = load('data_JT/migration_baseline.mat');
migration_has_reform = load('data_JT/migration_post_reform.mat');

% Unique cities
cityList = intersect(tradeMat.citygb(:),migration_cityid.CityId(:));
N = length(cityList)+1; % Accounting for RoW

% Look up index
[~,cityInTrade] = ismember(cityList,tradeMat.citygb);
cityInTrade = [cityInTrade;length(tradeMat.citygb)+1];    % Attach RoW
[~,cityInMigration] = ismember(cityList,migration_cityid.CityId);

% Extract trade and migration flows; calculate respective shares
X_od = tradeMat.X_od(cityInTrade,cityInTrade);
L_od = migration_has_reform.City_UtoU_s(cityInMigration,cityInMigration);
% Following the convention of CDP (2019), kappa_ni is trade cost from i to
% n, so take transpose
dot_kappa_ni = ( tradeMat.tau_od_2010(cityInTrade,cityInTrade) ./ tradeMat.tau_od_2000(cityInTrade,cityInTrade) )';
% Following the convention of CDP (2019), tau_ni is migration cost from n to i
dot_tau_ni = ( migration_has_reform.MigrationCost_UtoU_s(cityInMigration,cityInMigration) ./ migration_no_reform.MigrationCost_UtoU_s(cityInMigration,cityInMigration) );
dot_tau_ni = [
    dot_tau_ni,ones(N-1,1);
    ones(1,N-1),1
    ];  % Accounting for RoW

% X_n_0 is the expenditure normalized to one
X_n_0 = sum(X_od,1);
pi_ni_0 = permute(X_od ./ X_n_0, [2,1]);
X_n_0 = X_n_0 / sum(X_n_0); % normalization

% Calcualte the migration share;
mu_ni_m1 = L_od ./ sum(L_od,2);
mu_ni_m1 = mu_ni_m1^(1/50);   % Jingting's data is life-time migration flow; convert to one-period flow
mu_ni_m1 = [
    mu_ni_m1, zeros(N-1,1);
    zeros(1,N-1), 1
    ];  % Account for RoW
% L_n_0 is the aggregate labor normalized to one
% Scale to accomodate RoW
%{
L_d_domestic = sum(L_od,1)';
L_RoW = sum(L_d_domestic) / sum(tradeMat.L_d(1:end-1)) * tradeMat.L_d(end);
L_n_0 = [L_d_domestic;L_RoW] / (sum(L_d_domestic)+L_RoW);
%}
L_n_0 = tradeMat.L_d / sum(tradeMat.L_d);

%%%%%%%%%%%%%%%%% Inspects costs and plot %%%%%%%%%%%
change_domestic_trade_costs = log(dot_kappa_ni(1:end-1,1:end-1));
change_intl_trade_costs = log(dot_kappa_ni(1:end-1,end));
figure; hold on;
histogram(change_domestic_trade_costs,'normalization','pdf');
histogram(change_intl_trade_costs,'normalization','pdf');
legend({'$\Delta$ log domestic trade costs','$\Delta$ log int''l trade costs'},'FontSize',13);
title('Change in bilateral trade costs');
ylabel('Densities');
print('figures/delta_trade_costs.png','-dpng');

change_migration_costs = log(dot_tau_ni);
figure; hold on;
histogram(change_migration_costs,'normalization','pdf');
legend({'$\Delta$ log migration costs'},'FontSize',13);
title('Change in bilateral migration costs');
ylabel('Densities');
print('figures/delta_migration_costs','-dpng');

%%%%%%%%% Transition period
T = 100;

%%%%%%%%%% Input-output tables
gamma_n = 0.5*ones(N,1);
gamma_tilde_n = 0.5*ones(N,1);

%%%%%%%%%%%%%% Shocks
% baseline
dot_A_n_t = ones(N,T);
dot_kappa_ni_t = ones(N,N,T);
dot_kappa_ni_t(:,:,2) = dot_kappa_ni;
dot_tau_ni_t = ones(N,N,T);
dot_tau_ni_t(:,:,2) = dot_tau_ni;

% counterfactuals
hat_A_n_t = ones(N,T);
hat_kappa_ni_t = ones(N,N,T);
hat_kappa_ni_t(:,:,2) = 1./dot_kappa_ni;
hat_tau_ni_t = ones(N,N,T);
hat_tau_ni_t(:,:,2) = 1./dot_tau_ni;


clear params;
data = v2struct;
end