%this file is to generate the artifical fundamentals and parameters
%for model 1-3

%%
%number of regions
N=60; 

%artifical fundamentals
rng('default') %random seed
L_i = rand(N,1); %exogenous labor endowment (for model w/o migration)
T_i = rand(N,1); %productivity
A_n = rand(N,1); %amenity
L_i0=sum(L_i)./N; %initial labor distribution (for model w/ migration)

%(asymmetric) artifical trade cost
d_ni= rand(N,N)+1;
d_ni=d_ni.*(1-eye(size(d_ni)))+eye(size(d_ni));

%(asymmetric) artifical migration cost
tau_ni= 5*rand(N,N)+5;
tau_ni=tau_ni.*(1-eye(size(d_ni)))+eye(size(tau_ni));

%migration els.
epsilon=3;

%trade els.
theta=4;

%substitution els.
Sigma=4;

%labor share in CD production function
Alpha=0.8;



data=v2struct;





