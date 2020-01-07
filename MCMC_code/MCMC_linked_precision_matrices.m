function [C_save, Sig_save, adj_save, Phi_save, ar_Phi] = ...
  MCMC_linked_precision_matrices(S, n, C, Phi, v0, v1, lambda, pii, ...
      burnin, nmc, disp)
% Generate posterior sample of precision matrices, graphs, and Phi matrix
% Note that code for graph sampling is based on that provided by Wang (2014)
% Input parameters:
%   S: p x p x K array of sample covariance matrix for each group
%   n: vector of sample sizes for each group
%   C: p x p x K array of initial precision matrices for each group
%   Phi: initial value for correlation matrix linking edge values across groups
%   v0: prior parameter for standard deviation for included edges (large)
%   v1: prior parameter for standard deviation for non-included edges (small)
%   lambda: prior parameter for exponential on diagonal precision matrix entries
%   pii: prior probability 
%   burnin: number of burn-in iterations
%   nmc: number of MCMC sample iterations saved
%   disp: T/F for whether to print progress to screen
% Output parameters:
%   C_save: p x p x K x nmc sample of precision matrix
%   Sig_save: p x p x K x nmc sample of covariance matrix
%   adj_save: p x p x K x nmc sample of adjacency matrix
%   Phi_save: K x K x nmc sample of cross-graph correlation matrix
%   ar_Phi: acceptance rate for Phi

% K is number of sample groups
K = size(Phi, 1);

% p is number of variables
[p] = size(S, 1);

% Set up matrices for return values
C_save = zeros(p, p, K, nmc);
Sig_save = C_save;
adj_save = C_save;
Phi_save = zeros(K, K, nmc);

% Acceptance rate for Phi
ar_Phi = 0;

% Initial adjacency matrices as defined by intial precision matrices
adj = abs(C) > 1e-5;

% Elementwise product zeros out elements of C that are close to 0
C = C .* adj;

% Initial value for covariance matrix
Sig = zeros(p, p, K);
for k = 1:K
    Sig(:, :, k) = inv(C(:, :, k));
end

% Setup parameters for graph/precision matrix sampling
V0 = v0 * ones(p);
V1 = v1 * ones(p);

% Initial value for tau
tau = repmat(V1, [1, 1, K]);

ind_noi_all = zeros(p-1, p);
for i = 1:p
    if i == 1
        ind_noi = [2:p]';
    elseif i == p
        ind_noi = [1:p-1]';
    else
        ind_noi = [1:i-1,i+1:p]';
    end
    ind_noi_all(:, i) = ind_noi; 
end

ind_nok_all = zeros(K-1, K);
for i = 1:K
    if i == 1
        ind_nok = [2:K]';
    elseif i == K
        ind_nok = [1:K-1]';
    else
        ind_nok = [1:i-1,i+1:K]';
    end
    ind_nok_all(:, i) = ind_nok; 
end

% Perform MCMC sampling
for iter = 1:burnin + nmc
    if (disp && mod(iter, 500) == 0)
        fprintf('iter = %d\n', iter);
    end

    % Update graph and precision matrix for each group
    for k = 1:K
	ind_nok = ind_nok_all(:, k);

        for i = 1:p
            ind_noi = ind_noi_all(:, i);

            % Equivalent to v_12 vector in Wang (2014) paper
            tau_temp = tau(ind_noi, i, k);

            % Compute m vector and diagonal entries in D matrix
            D_diag = repmat(NaN, p-1, 1);
            m = repmat(NaN, p-1, 1);
            for cur_ind = 1:p-1
                cur_ind_noi = ind_noi(cur_ind);

                % Note that stored value of tau is in fact squared version
                % Use vector of nu_ij values across K groups as diagonal elements of matrix
                diag_nuij = diag(sqrt(squeeze(tau(cur_ind_noi, i, :))));

                % Cross graph relationship matrix for current edge i.e. Theta_ij in paper
                theta_xp = diag_nuij * Phi * diag_nuij;
                theta_xp11 = theta_xp(ind_nok, ind_nok);
                theta_xp12 = theta_xp(ind_nok, k);
                prod_xp = theta_xp12' / theta_xp11;
                m(cur_ind) = prod_xp * squeeze(C(cur_ind_noi, i, ind_nok));
                D_diag(cur_ind) = tau_temp(cur_ind) - prod_xp * theta_xp12;
            end 
        
            Sig11 = Sig(ind_noi, ind_noi, k);
            Sig12 = Sig(ind_noi, i, k);
        
            invC11 = Sig11 - Sig12 * Sig12' / Sig(i, i, k);
        
            % Compute core of C matrix (before inversion)
            Ci = (S(i, i, k) + lambda) * invC11 + diag(1 ./ D_diag);
            
            % Avg with its transpose to ensure symmetry
            Ci = (Ci + Ci') ./ 2;

            % Compute Cholesky decomposition of Ci
            Ci_chol = chol(Ci);

            % Compute a vector
            a = diag(1 ./ D_diag) * m - S(ind_noi, i, k);

            % Use this to compute Ca where C = Ci^-1
            mu_i = Ci_chol \ (Ci_chol' \ a);

            % Generate vector from multivariate normal distribution with mean mu_i and
            % covariance matrix Ci
            % Beta here is equivalent to u vector in Wang paper
            beta = mu_i + Ci_chol \ randn(p-1, 1);
         
            % Update off-diagonal elements for current row/col
            C(ind_noi, i, k) = beta;
            C(i, ind_noi, k) = beta;
        
            % Generate sample from Gamma distribution with parameters a_gam and b_gam
            % Gam here is equivalent to v in Wang paper
            a_gam = 0.5 * n(k) + 1;
            b_gam = (S(i, i, k) + lambda) * 0.5;
            gam = gamrnd(a_gam, 1 / b_gam);
            
            % Sampled value was after change of variables. Transform to original var
            % of diagonal entry for current row/col
            c = beta' * invC11 * beta;
            C(i, i, k) = gam + c;
        
            % Updating Covariance matrix according to one-column change of precision matrix
            invC11beta = invC11 * beta;
            Sig(ind_noi, ind_noi, k) = invC11 + invC11beta * invC11beta' / gam;
            Sig12 = -invC11beta / gam;
            Sig(ind_noi, i, k) = Sig12;
            Sig(i, ind_noi, k) = Sig12';
            Sig(i, i, k) = 1 / gam;  
        
            v0 = V0(ind_noi, i);
            v1 = V1(ind_noi, i);
          
            % Bernoulli update to adjacency matrix
            w1 = -0.5 * log(v0) - 0.5 * beta .^ 2 ./ v0 + log(1 - pii);
            w2 = -0.5 * log(v1) - 0.5 * beta .^ 2 ./ v1 + log(pii);
        
            w_max = max([w1,w2], [], 2);
        
            w = exp(w2-w_max) ./ sum(exp([w1,w2] - repmat(w_max,1,2)), 2);
        
            z = (rand(p-1, 1) < w);
        
            v = v0;
            v(z) = v1(z);
        
            tau(ind_noi, i, k) = v;
            tau(i, ind_noi, k) = v;
        
            adj(ind_noi, i, k) = z;
            adj(i, ind_noi, k) = z;  
        end
    end
    
    % Udpate cross-graph correlation matrix Phi
    
    % Compute Inverse-Wishart degrees of freedom
    IW_df = p * (p - 1) / 2;
    
    % Compute Inverse-Wishart scale matrix
    diag_V = repmat(0, K, 1);
    for k = 1:K
        diag_V(k) = 1 ./ sqrt(sum(sum(tril(C(:, :, k) .^ 2))));
    end
    V = diag(diag_V);
    
    IW_scale = zeros(K, K);
    for i = 1:p-1
        for j = i+1:p
            omega_ij = squeeze(C(i, j, :));
            e_ij = V * omega_ij;
            diag_nuij_inv = diag(1 ./ sqrt(squeeze(tau(i, j, :))));
            IW_scale = IW_scale + diag_nuij_inv * e_ij * e_ij' * diag_nuij_inv;
        end
    end
    
    % Sample covariance matrix
    % This gives us updated values for both Phi and V
    Psi_prop = iwishrnd(IW_scale, IW_df);
    V_inv_prop = diag(1 ./ sqrt(diag(Psi_prop)));
    Phi_prop = V_inv_prop * Psi_prop * V_inv_prop;
    
    % Compute MH acceptance ratio on log scale
    log_mh_ratio_Phi = (K + 1) / 2 * (logdet(Phi_prop) - logdet(Phi));
    
    % Accept with probability r
    if (log(rand(1)) < log_mh_ratio_Phi) 
        Phi = Phi_prop;
        ar_Phi = ar_Phi + 1 / (burnin + nmc);
    end

    % Retain values for posterior sample
    if iter > burnin
        C_save(:, :, :, iter-burnin) = C(:, :, :);
        Sig_save(:, :, :, iter-burnin) = Sig(:, :, :);
        adj_save(:, :, :, iter-burnin) = adj(:, :, :);
        Phi_save(:, :, iter-burnin) = Phi(:, :);
    end
end

