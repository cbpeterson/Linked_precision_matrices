clear;
%%%%%
%Add path to input files and MCMC code
addpath './MCMC_code';
addpath './Helper_functions';
addpath './Case_study/Code';

% Four graphs with simulated data
K = 4;
nHC_hp = 143;
nHC_Mem = 145;
nMCI = 148;
nAD = 148;

% Fix seed for random number generator
rng(3487841);

p = 100;
K = 4;

% (h in code) =  (h in paper )^2
h = 10^2;
    
% (v0 in code) = (v0 in paper)^2
v0 = 0.5;
    
% (v1 in code) = (v1 in paper)^2
v1 = 15;
    
lambda = 1;
pii = 2 / (p - 1);
    
% Joint inference code expects a vector of length K
pii = [pii pii pii pii];

% Prior parameters for gamma slab of mixture prior
alpha = 2;
beta = 5;
    
% Parameters for prior on nu which affects graph sparsity
a = 1;
b = 4;
    
% Parameter for Bernoulli prior on indicator of graph relatedness
w = .5;  

% Number of burn-in and post burn-in iterasations
burnin  = 10000;
nmc = 10000;

% Current iteration
% 
% fprintf('iter = %d\n', iter);
    
% Set up input S for each group with sample size n
n1 = nHC_hp;
cur_filename = strcat('./Case_study/Input/fake_hpHC.csv');
Y1 = readtable(cur_filename);
Y1 = table2array(Y1);
%%%%%%%
S1 = Y1' * Y1;

n2 = nHC_Mem;
cur_filename = strcat('./Case_study/Input/fake_HC.csv');
Y2 = readtable(cur_filename);
Y2 = table2array(Y2);
%%%%%%%
S2 = Y2' * Y2;

    
n3 = nMCI;
cur_filename = strcat('./Case_study/Input/fake_MCI.csv');
Y3 = readtable(cur_filename);
Y3 = table2array(Y3);
%%%%%%%
S3 = Y3' * Y3;
                                                                                                                                                                                                                                                    
n4 = nAD;
cur_filename = strcat('./Case_study/Input/fake_AD.csv');
Y4 = readtable(cur_filename);
Y4 = table2array(Y4);
S4 = Y4' * Y4;


    
% Initial value for precision matrix and covariance matrices
C1 = eye(p);
C2 = eye(p);
C3 = eye(p);
C4 = eye(p);
Sig1 = eye(p);
Sig2 = eye(p);
Sig3 = eye(p);
Sig4 = eye(p);

% Initial values for Theta and nu
Theta = zeros(K);
nu = zeros(p, p) - 1;

% Graph Specific Hyperparameters
V0_1 = v0 * ones(p);
V1_1 = v1 * ones(p);
V0 = cat(3, V0_1, V0_1, V0_1, V0_1);
V1 = cat(3, V1_1, V1_1, V1_1, V1_1);

% Call MCMC sampler
Sig = cat(3, Sig1, Sig2, Sig3, Sig4);
S = cat(3, S1, S2, S3, S4);
C = cat(3, C1, C2, C3, C4);
disp = true;
[C_save, Sig_save, adj_save, Theta_save, ar_gamma, ar_theta, ...
    nu_save, ar_nu] = ...
    MCMC_multiple_graphs_SSVS_Final(Theta, Sig, V0, V1, lambda, pii, ...
    [n1, n2, n3, n4], S, C, nu, alpha, beta, a, b, w, ...
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
[C_save_resample, Sig_save_resample] = ...
    MCMC_resample_C_for_fixed_graph(median_model, S, [n1, n2, n3, n4], C, Phi, ...
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

%% write files to be used in later analysis %%
csvwrite('./Case_study/Output/hpHC_joint.csv',ppi(:,:,1))
csvwrite('./Case_study/Output/HC_joint.csv',ppi(:,:,2))
csvwrite('./Case_study/Output/MCI_joint.csv',ppi(:,:,3))
csvwrite('./Case_study/Output/AD_joint.csv',ppi(:,:,4))



