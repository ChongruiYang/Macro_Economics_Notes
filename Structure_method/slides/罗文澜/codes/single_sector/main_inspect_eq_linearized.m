% This script inspects the transition path
clear;
linearized = load('mat/linearized.mat');

northeastProvList = [
    23  % Heilongjiang
    22  % Jilin
    21  % Liaoning
    ];

%%%%%%%%%%%%%%%%% Construct model statistics %%%%%%%%%%%%%%%%%%%%
isCityNortheast = ismember(floor(linearized.data.cityList/100),northeastProvList);
eqNameList = {'baseline','cf'};
legendList = {'Baseline','No Hukou Reform and No Expressway'};
lineSpecList = {'-','--',':'};
stats = struct;
N = linearized.data.N;
for i=1:length(eqNameList)
    eqName = eqNameList{i};
    stats.(eqName).ln_L_northeast_t = sum( ...
        linearized.(eqName).ln_L_n_t(isCityNortheast,:).*linearized.data.L_n_0(isCityNortheast)...
        ) ./ sum(linearized.data.L_n_0(isCityNortheast));
end
delta_ln_L_northeast_t = stats.baseline.ln_L_northeast_t - stats.cf.ln_L_northeast_t;

%%%%%%%%%%%%%%%% Population %%%%%%%%%%%%%%%
figure; hold on;
plot(2000:2099,delta_ln_L_northeast_t,'LineWidth',1.5);
xlim([2000,2030]);
legend(['$\Delta$ Log (population)'],'FontSize',11,'Location','Best','interpreter','latex');
xlabel('Year','interpreter','latex');
ylabel('$\Delta$ Log (population)','interpreter','latex');
title('Impact of Hukou and Expressway, from baseline to cf','interpreter','latex','fontsize',12);
print('figures/ln_L_northeast_linear.png','-dpng');


figure; hold on;
plot(2000:2099,exp(stats.baseline.ln_L_northeast_t),'LineWidth',1.5);
xlim([2000,2030]);
legend(['Baseline'],'FontSize',11,'Location','Best','interpreter','latex');
xlabel('Year','interpreter','latex');
title('Population relative to 2000','interpreter','latex','fontsize',12);
