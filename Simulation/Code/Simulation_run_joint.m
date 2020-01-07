% Add path to input files and MCMC code
addpath './Simulation/Input_files';
addpath './MCMC_code';
addpath './Helper_functions';

% Three graphs with simulated data
K = 3;
n = 150;
    
% Read in csv with precision matrices for each group
A1 = csvread('A1_sim_p100.csv');
A2 = csvread('A2_sim_p100.csv');
A3 = csvread('A3_sim_p100.csv');

% Number of variables = 100
p = size(A1, 1);
    
% Binary matrix indicating which edges should be present in G i.e. true
% adjancencies
adjTrue1 = abs(A1) > 0.001;
adjTrue2 = abs(A2) > 0.001;
adjTrue3 = abs(A3) > 0.001;

% Prior parameters for graph structure learning
% Note that the parameterization here is slightly different from that used in Wang (2014)

    
% (v0 in code) = (v0 in paper)^2
v0 = 0.05^2;
    
% (v1 in code) = (v1 in paper)^2
v1 = 0.5^2;
    
lambda = 1;
pii = 2 / (p - 1);
    
% Joint inference code expects a vector of length K
pii = [pii pii pii];

% Prior parameters for gamma slab of mixture prior
alpha = 2;
beta = 5;
    
% Parameters for prior on nu which affects graph sparsity
a = 1;
b = 16;
    
% Parameter for Bernoulli prior on indicator of graph relatedness
w = .5;  
    
% Number of burn-in and post burn-in iterations
burnin  = 10000;
nmc = 20000;

for run = 1:25
    fprintf('%d\n', run)
    datestr(now)
  
    % Fix seed for random number generator
    rng(3487841 + run);
    
    % Current iteration
    iter = run;
    
    fprintf('iter = %d\n', iter);
    
    % Set up input S for each group with sample size n
    n1 = n;
    cur_filename = strcat('Y1_iter', num2str(iter), '_n', num2str(n), '_lpm_p100.csv');
    Y1 = csvread(cur_filename);
    S1 = Y1' * Y1;
    
    n2 = n;
    cur_filename = strcat('Y2_iter', num2str(iter), '_n', num2str(n), '_lpm_p100.csv');
    Y2 = csvread(cur_filename);
    S2 = Y2' * Y2;
    
    n3 = n;
    cur_filename = strcat('Y3_iter', num2str(iter), '_n', num2str(n), '_lpm_p100.csv');
    Y3 = csvread(cur_filename);
    S3 = Y3' * Y3;
    
    % Initial value for precision matrix and covariance matrices
    C1 = eye(p);
    C2 = eye(p);
    C3 = eye(p);
    Sig1 = eye(p);
    Sig2 = eye(p);
    Sig3 = eye(p);
    
    % Initial values for Theta and nu
    Theta = zeros(K);
    nu = zeros(p, p) - 1;
    
    % Graph Specific Hyperparameters
    V0_1 = v0 * ones(p);
    V1_1 = v1 * ones(p);
    V0 = cat(3, V0_1, V0_1, V0_1);
    V1 = cat(3, V1_1, V1_1, V1_1);

    % Call MCMC sampler
    Sig = cat(3, Sig1, Sig2, Sig3);
    S = cat(3, S1, S2, S3);
    C = cat(3, C1, C2, C3);
    disp = true;
    [C_save, ~, adj_save, Theta_save, ar_gamma, ar_theta, ...
        nu_save, ar_nu] = ...
        MCMC_multiple_graphs_SSVS_Final(Theta, Sig, V0, V1, lambda, pii, ...
        [n1, n2, n3], S, C, nu, alpha, beta, a, b, w, ...
        burnin, nmc, disp);
    
    % PPIs for Theta (graph similarity measure)
    ppi_theta = mean(Theta_save ~= 0, 3);
    mean_theta = mean(Theta_save, 3);
    
    % Compute posterior probabilities of edge inclusion
    ppi = mean(adj_save, 4);
    
    % Take median model as selected graph for resampling conditional on graph
    median_model = ppi > 0.5;
    
    % Resample precision matrices conditional on graph and Phi
    % Fix Phi to identity matrix so that precision matrices are sampled
    % independently
    Phi = eye(K);
    [C_save_resample, ~] = ...
        MCMC_resample_C_for_fixed_graph(median_model, S, [n1, n2, n3], C, Phi, ...
        v0, v1, lambda, burnin, nmc, disp);
    
    % Estimate C (concentration/precision matrix) as MCMC avg
    est_C = mean(C_save, 4);
    est_C_resample = mean(C_save_resample, 4);
    
    % Calculate ppi of difference based on full MCMC sample
    ppi_diff = zeros(p * K, p * K);
    for k1 = 1:(K - 1)
        for k2 = (k1 + 1):K
            cur_adj_diff = abs(squeeze(adj_save(:, :, k1, :) - adj_save(:, :, k2, :)));
            cur_diff_ppi = mean(cur_adj_diff, 3);
            
            start_row = (k1 - 1) * p + 1;
            end_row = k1 * p;
            start_col = (k2 - 1) * p + 1;
            end_col = k2 * p;
            
            ppi_diff(start_row:end_row, start_col:end_col) = cur_diff_ppi;
        end
    end
    
    % Compute performance of edge selection and precision matrix estimation
    [tpr_cur, fpr_cur, mcc_cur, auc_cur, tpr_diff, fpr_diff, mcc_diff, auc_diff, ...
        fl] = get_perf(...
        cat(3, adjTrue1, adjTrue2, adjTrue3), ...
        ppi, ppi_diff, ...
        cat(3, A1, A2, A3), ...
        est_C_resample);
    
    % Save performance summary and edge PPIs
    cur_file = strcat('perf_iter', num2str(iter), '_n', num2str(n), '_joint_sim_p100.csv');
    csvwrite(cur_file, {tpr_cur, fpr_cur, mcc_cur, auc_cur, tpr_diff, fpr_diff, ...
        mcc_diff, auc_diff, fl});
    cur_file = strcat('PPI_iter', num2str(iter), '_n', num2str(n), '_joint_sim_p100.csv');
    csvwrite(cur_file, ppi);
    datestr(now)
end
