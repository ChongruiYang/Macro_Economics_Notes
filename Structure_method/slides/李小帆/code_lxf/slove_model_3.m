%this function is to solve model-3
%multiple regios
%one sector with  intermediate input
%with migration
%two step:
% 1).given labor distribution after migration, solve trade block (using solve_model_2.m)
% 2).combine initial labor dist. and real income from step-1, update guess of labor distribution
% 3).return to step-1 until converge
%%
function eqlm=slove_model_3(data);

v2struct(data);

%initial guess of L_i (labor distribution after migration)
L_i=L_i0;

%%%%%%
%iteration on labor dist. L_i
%%%%%%
diff3=1;
while diff3>1e-4;
    
    %%%solve trade block given guess of labor distribution
    data.L_i=L_i; %labor dist. for solving trade block
    trade_eqlm=slove_model_2(data); %solve trade block
    v2struct(trade_eqlm);
    
    U_n=trade_eqlm.W; %real income in region n
    
    %%%calculate new labor distribution combine initial guess and realincome
    %migration share
    lambda_ni=A_n.*(tau_ni.^(-1).*U_n).^(epsilon)./...
        sum(A_n.*(tau_ni.^(-1).*U_n).^(epsilon));
    %sum(lambda_ni);
    
    %implied labor distribution
    L_n_new=sum(lambda_ni.*L_i0',2);
    L_i_new=L_n_new;
    
    diff3=max(abs(L_i-L_i_new));  %converge condition
    L_i=0.5*L_i+0.5*L_i_new; %%%update guess of L_i until converge
end



eqlm=v2struct;






end