%this file is to generate parameters and solve equilibrium 

%path-set
addpath(genpath(pwd));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc
gen_parameters; 
eqlm=slove_model_1(data);
v2struct(eqlm);

%%%check trade balance
exp_i=sum((1-eye(N)).*pi_ni.*(w_n.*L_n))';
imp_n=sum((1-eye(N)).*pi_ni.*(w_n.*L_n),2);
imp_i=imp_n;
scatter(exp_i,imp_i)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check ACR(2012) sufficient statistics. of gain from trade
clear;clc
gen_parameters;
data.d_ni(d_ni>1)=Inf;
eqlm1=slove_model_1(data);

data.d_ni=ones(N,N);
eqlm2=slove_model_1(data);

w_hat=eqlm2.W./eqlm1.W;
openness=diag(eqlm2.pi_ni).^(-1/theta);
scatter(openness,w_hat)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc
gen_parameters;
eqlm=slove_model_2(data);

v2struct(eqlm);
X_i=w_i.*L_i+((1-Alpha)/Alpha).*w_i.*L_i;
X_n=X_i;

exp_i=sum((1-eye(N)).*pi_ni.*X_n)';
imp_n=sum((1-eye(N)).*pi_ni.*X_n,2);
imp_i=imp_n;

scatter(exp_i,imp_i)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc
gen_parameters;
%data.tau_ni(data.tau_ni>1)=Inf;
eqlm=slove_model_3(data);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model-4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc
gen_parameters2; data.d_nij(data.d_nij>1)=1;
eqlm=slove_model_4(data);

v2struct(eqlm);
exp_ij=permute(sum((1-eye(N)).*pi_nij.*reshape(X_nj,[N,1,J])),[2,3,1]);
imp_nj=permute(sum((1-eye(N)).*pi_nij.*reshape(X_nj,[N,1,J]),2),[1,3,2]);
imp_ij=imp_nj;

EXP_i=sum(exp_ij,2);
IMP_i=sum(imp_ij,2);
scatter(EXP_i,IMP_i)

