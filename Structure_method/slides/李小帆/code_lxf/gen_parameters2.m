%this file is to generate the artifical fundamentals and parameters
%for model 4

%%
%number of regions
N=60; 
%number of sectors
J=20;

%artifical fundamentals
rng('default') %random seed
L_i = rand(N,1); %exogenous labor endowment in region i (for model w/o migration)
T_ij = rand(N,J);  %productivity of sector j in region i

%(asymmetric) sector-level trade cost
d_nij= rand(N,N,J)+1;
d_nij=d_nij.*(1-eye(N))+eye(size(N));

%artifical Input-Output linkage
gamma_nj=0.4+(0.6)*rand(N,J); %labor share of sector j in region n (above 0.4)
gamma_nkj=rand(N,J,J); %share of input from k to j in region n
%normalize so as to labor share + input share = 1
gamma_nkj=(gamma_nkj./sum(gamma_nkj,2)).*(1-reshape(gamma_nj,[N,1,J]));
%test: permute(sum(gamma_nkj,2),[1,3,2])+gamma_nj;

%consumption share
Alpha_nj=0.4+(0.6)*rand(N,J);
Alpha_nj=Alpha_nj./sum(Alpha_nj,2);

%trade els.
theta_j=4*ones(J,1);

%substitution els.
Sigma_j=4*ones(J,1);



data=v2struct;





