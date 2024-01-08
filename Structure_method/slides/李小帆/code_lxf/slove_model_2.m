%this function is to solve model-2
%multiple regios
%one sector with  intermediate input
%no migration

%%
function eqlm=slove_model_2(data);

v2struct(data);

%initial guess of w_i and P_i
w_i=ones(N,1);
P_i=ones(N,1);

%constant term in price index
Ga=gamma((theta+1-Sigma)/theta)^(1/(1-Sigma));

%%%%%%
%iteration on wage w_i
%%%%%%
diff2=1;
while diff2>1e-4
    
    %%%%%%
    %iteration on price P_i
    %%%%%%
    diff1=1;
    while diff1>1e-4
        %%%calculate price (given wage and price guess)
        c_i=w_i.^Alpha.*P_i.^(1-Alpha);  %unit cost
        P_n=Ga*sum(T_i'.*(c_i'.*d_ni).^(-theta),2).^(-1/theta); %implied price
        P_i_new=P_n;
        
        diff1=max(abs(P_i-P_i_new)); %converge condition
        P_i=P_i_new; %update guess of price until converge
    end
    
    %trade share
    pi_ni=T_i'.*(c_i'.*d_ni).^(-theta)./sum(T_i'.*(c_i'.*d_ni).^(-theta),2);
    
    %%%calculate expenditure and income (given wage)
    w_n=w_i;
    L_n=L_i;
    expend_i=w_i.*L_i; %expenditure
    income_i= sum(pi_ni.*(w_n.*L_n))';  %income
    
    diff2=max(abs(expend_i-income_i));  %converge condition
    
    %%%update guess of w_i until converge
    w_i=0.5*w_i+0.5*w_i.*(income_i./expend_i).^(0.5);
    w_i=w_i/w_i(1); %normalization
end

P_i=P_n;
W=w_i./P_i; %real income or welfare



eqlm=v2struct;
end


