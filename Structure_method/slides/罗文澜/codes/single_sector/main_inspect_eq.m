% This script inspects the transition path
% clear;
eqList.baseline = load('mat/baseline.mat');
eqList.cf_no_hukou_no_expressway = load('mat/cf.mat');
eqList.cf_no_hukou = load('mat/cf_no_hukou.mat');

northeastProvList = [
    23  % Heilongjiang
    22  % Jilin
    21  % Liaoning
    ];

%%%%%%%%%%%%%%%%% Construct model statistics %%%%%%%%%%%%%%%%%%%%
isCityNortheast = ismember(floor(eqList.baseline.data.cityList/100),northeastProvList);
eqNameList = {'baseline','cf_no_hukou','cf_no_hukou_no_expressway'};
legendList = {'Baseline','No Hukou Reform','No Hukou Reform and No Expressway'};
lineSpecList = {'-','--',':'};
stats = struct;
N = eqList.baseline.data.N;
for i=1:length(eqNameList)
    eqName = eqNameList{i};
    if i==1
        stats.(eqName).L_northeast_t = sum(eqList.(eqName).eq.L_n_t(isCityNortheast,:));
        stats.(eqName).omega_n_t = cumprod(eqList.(eqName).eq.dot_omega_n_t,2);
        stats.(eqName).omega_northeast_t = sum(stats.(eqName).omega_n_t(isCityNortheast,:) .* eqList.(eqName).eq.L_n_t(isCityNortheast,:)) ./ stats.(eqName).L_northeast_t;
        stats.(eqName).omega_national_t = sum(stats.(eqName).omega_n_t(1:N-1,:) .* eqList.(eqName).eq.L_n_t(1:N-1,:)) ./ sum(eqList.(eqName).eq.L_n_t(1:N-1,:));
        
    else
        stats.(eqName).L_northeast_t = sum(eqList.(eqName).eq.prime_L_n_t(isCityNortheast,:));
        stats.(eqName).omega_n_t = cumprod(eqList.baseline.eq.dot_omega_n_t .* eqList.(eqName).eq.hat_omega_n_t,2);
        stats.(eqName).omega_northeast_t = sum(stats.(eqName).omega_n_t(isCityNortheast,:) .* eqList.(eqName).eq.prime_L_n_t(isCityNortheast,:)) ./ stats.(eqName).L_northeast_t;
        stats.(eqName).omega_national_t = sum(stats.(eqName).omega_n_t(1:N-1,:) .* eqList.(eqName).eq.prime_L_n_t(1:N-1,:)) ./ sum(eqList.(eqName).eq.prime_L_n_t(1:N-1,:));
    end
    stats.(eqName).omega_rel_t = stats.(eqName).omega_northeast_t ./ stats.(eqName).omega_national_t;
end

%%%%%%%%%%%%%%%%%% Read data statistics %%%%%%%%%%
dataStats = load('data_Luo/northeast_data.mat');


%%%%%%%%%%%%%%%% Population %%%%%%%%%%%%%%%
figure; hold on;
plot(dataStats.yearList(2:end),dataStats.rel_pop(2:end) / dataStats.rel_pop(2),'k','LineWidth',1.5);
for i=1:length(eqNameList)
    eqName = eqNameList{i};
    plot(2000:2099,stats.(eqName).L_northeast_t/ stats.baseline.L_northeast_t(1),lineSpecList{i} ,'LineWidth',1.5);
end
xlim([2000,2030]);
legend(['Data',legendList],'FontSize',11,'Location','Best','interpreter','latex');
xlabel('Year','interpreter','latex');
ylabel('Relative to 2000','interpreter','latex');
title('Ratio of Northeast Population to National Population (Year 2000=1)','interpreter','latex','fontsize',12);
print('figures/rel_population_ratio.png','-dpng');

%%%%%%%%%%%%%% Log change %%%%%%%%%%
figure; hold on;
plot(2000:2099,log(stats.baseline.L_northeast_t) - log(stats.cf_no_hukou_no_expressway.L_northeast_t),'LineWidth',1.5);
xlim([2000,2030]);
legend(['$\Delta$ Log (population)'],'FontSize',11,'Location','Best','interpreter','latex');
xlabel('Year','interpreter','latex');
ylabel('$\Delta$ Log (population)','interpreter','latex');
title('Impact of Hukou and Expressway, from baseline to cf','interpreter','latex','fontsize',12);
print('figures/ln_L_northeast_nonlinear.png','-dpng');

%%%%%%%%%%% GDP per capita %%%%%%%%%%%%%%
figure; hold on;
plot(dataStats.yearList(2:end),dataStats.rel_real_gdp(2:end) / dataStats.rel_real_gdp(2),'k','LineWidth',1.5);
for i=1:length(eqNameList)
    eqName = eqNameList{i};
    plot(2000:2099,stats.(eqName).omega_rel_t/ stats.baseline.omega_rel_t(1),lineSpecList{i} ,'LineWidth',1.5);
end
xlim([2000,2030]);
legend(['Data',legendList],'FontSize',11,'Location','Best','interpreter','latex');
xlabel('Year','interpreter','latex');
ylabel('Relative to 2000','interpreter','latex');
title('Ratio of Northeast Real GDP p.c. to National (Year 2000=1)','interpreter','latex','fontsize',12);
print('figures/rel_gdp_pc_ratio.png','-dpng');
