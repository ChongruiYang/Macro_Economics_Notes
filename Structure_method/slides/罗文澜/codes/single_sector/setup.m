function params = setup
theta = 4;
beta = 0.5;
nu = 3*beta;

% algorithm
DAMP_P = 0.3;
DAMP_U = 0.5;
TOL_TEMPORARY_EQ = 1e-6;
TOL_DYNAMIC_EQ = 1e-5;
MINOR_PRINT_FREQ = inf;
PRINT_FREQ = 1;

params = v2struct;
end