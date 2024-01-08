%this function is to solve model-1
%multiple regions
%one sector with no intermediate input
%no migration

%%
function eqlm=slove_model_1(data);

v2struct(data);

%initial guess of w_i
w_i=ones(N,1);

%%%%%%
%iteration on wage w_i
%%%%%%
diff=1;
while diff>1e-4 %stop criterion
    
    %%%calculate trade share (given wage)
    c_i=w_i; %unit cost
    pinom_ni=T_i'.*(c_i'.*d_ni).^(-theta); %nominator of trade share pi_ni
    pidenom_n=sum(T_i'.*(c_i'.*d_ni).^(-theta),2); %denominator of trade share pi_ni
    pi_ni=pinom_ni./pidenom_n; %trade share
    %test: sum(pi_ni,2);
    
    %%%calculate expenditure and income (given wage)
    w_n=w_i;
    L_n=L_i;
    expend_i=w_i.*L_i; %expenditure
    income_i= sum(pi_ni.*(w_n.*L_n))'; %income
    
    diff=max(abs(expend_i-income_i));  %converge condition
    
    %%%update guess of w_i until converge
    w_i=0.5*w_i+0.5*w_i.*(income_i./expend_i).^(0.5);
    w_i=w_i/w_i(1); %normalization
end

Ga=gamma((theta+1-Sigma)/theta)^(1/(1-Sigma));
P_n=Ga*sum(T_i'.*(c_i'.*d_ni).^(-theta),2).^(-1/theta); %price index
P_i=P_n;
W=w_i./P_i; %real income or welfare



eqlm=v2struct;
end


