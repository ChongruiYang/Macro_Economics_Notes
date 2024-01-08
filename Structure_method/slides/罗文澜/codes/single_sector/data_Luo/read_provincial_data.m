% This script reads the provincial data
clear;

%%%%%%%%%%%%%% Population %%%%%%%%%%%%%%%%
startRow = 34;
endRow = startRow + 31;
sheetRead = readtable('provincial_data.xls','Range',['B',num2str(startRow),':W',num2str(endRow)]);
regionList = cellstr(sheetRead{:,1});
yearList = 1999:2019;
northeastList = {
    '辽宁省'
    '吉林省'
    '黑龙江省'
    };
dataList = sheetRead{:,2:end};
pop_national = dataList(1,:);
pop_northeast = sum(dataList(ismember(regionList,northeastList),:));
rel_pop_growth = pop_northeast(:,2:end) ./ pop_northeast(:,1:end-1) ...
    ./(pop_national(:,2:end) ./ pop_national(:,1:end-1));
rel_pop = pop_northeast ./ pop_national;

figure;
plot(yearList,rel_pop*100,'LineWidth',1.5);
xlim([1998,2020]);
xticks([2000:5:2020]);
% ylim([0.9,1.01]);
xlabel('Year');
ylabel('%');
title('Population of Northeast Region as Share of National Population','interpreter','latex');
print('rel_pop_northeast.png','-dpng');

%%%%%%%%%%%%%%%%%% Price level 
startRow = 98;
endRow = startRow + 31;
sheetRead = readtable('provincial_data.xls','Range',['B',num2str(startRow),':W',num2str(endRow)]);
regionList = cellstr(sheetRead{:,1});
yearList = 1999:2019;
northeastList = {
    '辽宁省'
    '吉林省'
    '黑龙江省'
    };
dataList = sheetRead{:,2:end};
cpi_national = dataList(1,:);
cpi_northeast = mean(dataList(ismember(regionList,northeastList),:));
price_national = cumprod([1,cpi_national(2:end)/100]);
price_northeast = cumprod([1,cpi_northeast(2:end)/100]);

%%%% GDP
sheetRead = readtable('gdp.xlsx','ReadVariableNames',false,'DatetimeType','text');
sheetSubset = sheetRead{51:71,2:end};
gdp_national = sum(sheetSubset,2)';
gdp_northeast = sum(sheetSubset(:,7:9),2)';

%%%%%%%%%%% real income
real_gdp_national = gdp_national ./ price_national ./ pop_national;
real_gdp_northeast = gdp_northeast ./ price_northeast ./ pop_northeast;
rel_real_gdp = real_gdp_northeast./real_gdp_national;
figure;
plot(yearList,real_gdp_northeast./real_gdp_national,'LineWidth',1.5);
xlim([1998,2020]);
xticks([2000:5:2020]);
xlabel('Year');
title('Average Real GDP Per Capita of Northeast Relative to National','interpreter','latex');
print('rel_real_gdp_northeast.png','-dpng');

save('northeast_data.mat','yearList','rel_pop','rel_real_gdp');