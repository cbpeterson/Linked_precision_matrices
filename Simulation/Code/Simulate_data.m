% Generate simulated data and write out as csv to allow comparison between
% methods implemented in R and those implemented in Matlab

% Set seed for random number generator
rng(41912);

% Add path to directory containing helper functions
addpath('./Simulation/Input_files');
addpath('./Helper_functions')

% Read in csv with precision matrices for each group created by script
% Set_up_precision_matrices.m
A1 = csvread('A1_sim_p100.csv');
A2 = csvread('A2_sim_p100.csv');
A3 = csvread('A3_sim_p100.csv');

% True covariance matrices are the inverse of these precision matrices
Cov1_True = inv(A1);
Cov2_True = inv(A2);
Cov3_True = inv(A3);

% Number of iterations per simulation (for each n setting)
niter = 25;

% Number of variables = 100
p = 100;

% Run for n = 150
n = 150;

% Sample sizes are equal in all three groups
n1 = n;
n2 = n;
n3 = n;

for iter = 1:niter    
    % Create Y for group 1 with sample size n1 and save as csv
    Y1 = rMNorm(zeros(p, 1), Cov1_True, n1)';
    cur_filename = strcat('Y1_iter', num2str(iter), '_n', num2str(n), '_lpm_p100.csv');
    csvwrite(cur_filename, Y1);
    
    % Create Y for group 2 with sample size n2 and save as csv
    Y2 = rMNorm(zeros(p, 1), Cov2_True, n2)';
    cur_filename = strcat('Y2_iter', num2str(iter), '_n', num2str(n), '_lpm_p100.csv');
    csvwrite(cur_filename, Y2);

    % Create Y for group 3 with sample size n3 and save as csv
    Y3 = rMNorm(zeros(p, 1), Cov3_True, n3)';
    cur_filename = strcat('Y3_iter', num2str(iter), '_n', num2str(n), '_lpm_p100.csv');
    csvwrite(cur_filename, Y3);
end
