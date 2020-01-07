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
v0 = 0.01^2;
    
% (v1 in code) = (v1 in paper)^2
v1 = 0.1^2;
    
lambda = 1;
pii = 2 / (p - 1);
    
% Number of burn-in and post burn-in iterations
burnin  = 10000;
nmc = 20000;

for run = 1:25
    fprintf('%d\n', run)
    datestr(now)

    % Fix seed for random number generator
    rng(3195700 + run);
    
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
    
    % Initial value for precision matrix
    C1 = eye(p);
    C2 = eye(p);
    C3 = eye(p);
    
    % Call MCMC sampler
    S = cat(3, S1, S2, S3);
    C = cat(3, C1, C2, C3);
    disp = true;
    [C_save, ~, adj_save] = ...
        MCMC_sep_precision_matrices(S, [n1, n2, n3], C, ...
        v0, v1, lambda, pii, burnin, nmc, disp);
    
    % Compute posterior probabilities of edge inclusion
    ppi = mean(adj_save, 4);
    
    % Take median model as selected graph for resampling conditional on graph
    median_model = ppi > 0.5;
    
    % Resample precision matrices conditional on graph and Phi (fixed to
    % identity matrix)
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
        fl_cur] = get_perf(...
        cat(3, adjTrue1, adjTrue2, adjTrue3), ...
        ppi, ppi_diff, ...
        cat(3, A1, A2, A3), ...
        est_C_resample);
    
    % Save performance summary, edge PPIs, and estimated Phi
    cur_file = strcat('perf_iter', num2str(iter), '_n', num2str(n), '_sep_sim_p100.csv');
    csvwrite(cur_file, {tpr_cur, fpr_cur, mcc_cur, auc_cur, tpr_diff, fpr_diff, ...
        mcc_diff, auc_diff, fl_cur});
    cur_file = strcat('PPI_iter', num2str(iter), '_n', num2str(n), '_sep_sim_p100.csv');
    csvwrite(cur_file, ppi);
    
    datestr(now)
end
