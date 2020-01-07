function [C_save, Sig_save, adj_save, Theta_save, ar_gamma, ar_theta, ...
    nu_save, ar_nu] = MCMC_multiple_graphs_SSVS_Final(Theta, Sig, V0, V1, lambda, pii, n, S, ...
    C, nu, alpha, beta, a, b, my_w, burnin, nmc, disp)
% Modified from code provided by Hao Wang to allow inference on multiple graphs
% Input parameters:
%   Theta: initial value for graph similarity matrix
%   Sig: initial guess for Sigma p x p x K matrix where C inverse = Sig
%   V0:  p x p x K matrix of the  small variance  components;  V0(i,j,n) is
%   the small variance of the precision element C(i,j) for group n
%   V1:  p x p x K matrix of the  large variance  components;  V1(i,j,n) is
%   the small variance of the precision element C(i,j) for group n
%   lambda: hyperparameter for the diagonal element Usually set to be 1;
%   pii: vector of prior marginal edge "inclusion probability" for each
%   group
%   n: vector of sample sizes for each group
%   S: p x p x K array of sample covariance matrix for each group
%   C: p x p x K array of initial precision matrices for each group
%   nu: p x p matrix of initial values for nu
%   alpha: shape parameter for Theta slab
%   beta: rate parameter for Theta slab
%   a: first parameter for prior on nu
%   b: second parameter for prior on nu
%   my_w: Bernoulli prior parameter for graph similarity indicators
%   burnin: number of burn-in iterations
%   nmc: number of MCMC sample iterations saved
%   disp: T/F for whether to print progress to screen
% Output parameters:
%   C_save: p x p x K x nmc sample of precision matrix
%   Sig_save: p x p x K x nmc sample of covariance matrix
%   adj_save: p x p x K x nmc sample of adjacency matrix
%   Theta_save: K x K x nmc sample of graph similarity matrix
%   ar_gamma: acceptance rate for between-model moves
%   ar_theta: acceptance rate for within-model moves
%   nu_save: p x p x nmc sample of edge-specific parameters
%   ar_nu: acceptance rate for nu


% K is number of sample groups
K = size(Theta, 1);

% p is number of variables
[p] = size(S, 1);

% Set up matrices for return values
C_save = zeros(p, p, K, nmc);
Sig_save = C_save;
adj_save = C_save;
Theta_save = zeros(K, K, nmc);
nu_save = zeros(p, p, nmc);

tau = V1;
 
 
 
ind_noi_all = zeros(p-1,p);
 
 
for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       ind_noi_all(:,i) = ind_noi;
       
end
 
 
pii_RB = zeros(p, p, K);
pii_mat = zeros(p, p, K);
for i = 1:K
    pii_mat(:,:,i) = pii(i);
end

%grid1 = [0:0.05:100];
%grid2 = [-60:0.05:60];
 
%%logd_offdiag_save = zeros(nmc,length(grid2));
%%logd_diag_save  = zeros(nmc,length(grid1));

% Acceptance rate for gamma based on number of between model moves
ar_gamma = zeros(K, K);

% Acceptance rate for theta based on number of within model moves
n_within_model = zeros(K, K);
ar_theta = zeros(K, K);
ar_nu = zeros(p, p);

% Initial adjacency matrices as defined by intial precision matrices
adj = abs(C) > 1e-5;

% Proposal parameters for MH steps
alpha_prop = 1;
beta_prop = 1;
a_prop = 2;
b_prop = 4;

% Elementwise product zeros out elements of C that are close to 0
C = C .* adj;

% Perform MCMC sampling
for iter = 1:burnin + nmc
    if (disp && mod(iter, 500) == 0)
        fprintf('iter = %d\n', iter);
    end

    % Update graph and precision matrix for each group using code from Hao Wang
    for cur_graph = 1:K

        % Sample Sig and C
        for i = 1:p
            ind_noi = ind_noi_all(:,i);
 
       
      tau_temp = tau(ind_noi,i,cur_graph);
       
      Sig11 = Sig(ind_noi,ind_noi,cur_graph); 
      Sig12 = Sig(ind_noi,i, cur_graph);
      
      invC11 = Sig11 - Sig12*Sig12'/Sig(i,i, cur_graph);
      
      Ci = (S(i,i,cur_graph)+lambda)*invC11+diag(1./tau_temp);  

  
      Ci = (Ci+Ci')./2;       
      Ci_chol = chol(Ci);    
      mu_i = -Ci_chol\(Ci_chol'\S(ind_noi,i,cur_graph));
        epsilon = mu_i+ Ci_chol\randn(p-1,1);
        
        
        C(ind_noi, i, cur_graph) = epsilon;
        C(i, ind_noi, cur_graph) = epsilon;
        
        a_gam = 0.5*n(:,cur_graph)+1;
        b_gam = (S(i,i,cur_graph)+lambda)*0.5;
        gam = gamrnd(a_gam,1/b_gam);
         
        c = epsilon'*invC11*epsilon;
        
        C(i,i, cur_graph) = gam+c;
    
        % note epsilon is beta in original Hao Wang code
        
        
        %% Below updating Covariance matrix according to one-column change of precision matrix
        invC11epsilon = invC11*epsilon;
        
        Sig(ind_noi,ind_noi,cur_graph) = invC11+invC11epsilon*invC11epsilon'/gam;
        Sig12 = -invC11epsilon/gam;
        Sig(ind_noi,i,cur_graph) = Sig12;
        Sig(i,ind_noi,cur_graph) = Sig12';
        Sig(i,i,cur_graph) = 1/gam;
 
        
        
              
         v0 = V0(ind_noi,i,cur_graph);
         v1 = V1(ind_noi,i,cur_graph);
         
         mrf_sum=0;
         for foo = 1:K
             mrf_sum=mrf_sum+Theta(cur_graph,foo)*pii_mat(ind_noi,i,foo);
         end
         
 % Applying the MRF prior to probability of inclusion
        pii_mat(ind_noi,i,cur_graph)=exp(pii_mat(ind_noi,i,cur_graph).*(nu(ind_noi,i)+2*mrf_sum))./(1+exp(nu(ind_noi,i)+2*mrf_sum));     
                
        w1 = -0.5*log(v0) -0.5*epsilon.^2./v0+log(1-pii_mat(ind_noi,i,cur_graph));
        w2 = -0.5*log(v1) -0.5*epsilon.^2./v1+log(pii_mat(ind_noi,i,cur_graph));
               
        %w1 = -0.5*log(v0) -0.5*epsilon.^2./v0+log(1-pii(:,cur_graph));
        %w2 = -0.5*log(v1) -0.5*epsilon.^2./v1+log(pii(:,cur_graph));
        
        w_max = max([w1,w2],[],2);
 
        w = exp(w2-w_max)./sum(exp([w1,w2]-repmat(w_max,1,2)),2);
 
        z = (rand(p-1,1)<w);
        
        
        v = v0;
        v(z) = v1(z);
        
        tau(ind_noi, i, cur_graph) = v;        
        tau(i, ind_noi, cur_graph) = v;
 
        adj(ind_noi, i, cur_graph) = z;
        adj(i, ind_noi, cur_graph) = z;
            
        end

        if iter > burnin
            pii_RB = pii_RB + pii_mat/nmc; 
            C_save(:, :, cur_graph, iter-burnin) = C(:, :, cur_graph);
            adj_save(:, :, cur_graph, iter-burnin) = adj(:, :, cur_graph);
            Sig_save(:,:,cur_graph, iter-burnin) = Sig(:,:,cur_graph); 
            
            
        %    logd_offdiag_save(iter-burnin,:) = logd_offdiag;
        %    logd_diag_save(iter-burnin,:) = logd_diag;
        
        end
    end

    
     %       logd_max = max(logd_offdiag_save);
     %        d_offdiag = exp(logd_max).*mean(exp(logd_offdiag_save-repmat(logd_max,nmc,1)));
 
 
     %       logd_max = max(logd_diag_save);
     %       d_diag = exp(logd_max).*mean(exp(logd_diag_save-repmat(logd_max,nmc,1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Update the parameters for network relatedness
    for k = 1:K-1
        for m = k+1:K
            % Between model move
            if Theta(k, m) == 0
                theta_prop = gamrnd(alpha_prop, beta_prop);
            else
                theta_prop = 0;
            end
            
            Theta_prop = Theta;
            Theta_prop(k, m) = theta_prop;
            Theta_prop(m, k) = theta_prop;
            
            % Get terms that are a sum over all edges on log scale
            sum_over_edges = 0;
            for i = 1:p-1
                for j = i+1:p
                    sum_over_edges = sum_over_edges + ...
                        log(calc_mrf_C(Theta, nu(i, j))) + 2 * ...
                        (theta_prop - Theta(m, k)) * adj(i, j, k) * adj(i, j, m) - ...
                        log(calc_mrf_C(Theta_prop, nu(i, j)));
                end
            end
            
            % Calculate MH ratio on log scale
            if theta_prop == 0
                log_ar = alpha_prop * log(beta_prop) - log(gamma(alpha_prop)) + ...
                    log(gamma(alpha)) - alpha * log(beta) - ...
                    (alpha - alpha_prop) * log(Theta(m, k)) + ...
                    (beta - beta_prop) * (Theta(m, k)) + sum_over_edges + ...
                    log(1 - my_w) - log(my_w);
            else
                log_ar = alpha * log(beta) - log(gamma(alpha)) + ...
                    log(gamma(alpha_prop)) - alpha_prop * log(beta_prop) - ...
                    (alpha - alpha_prop) * log(theta_prop) - ...
                    (beta - beta_prop) * theta_prop + sum_over_edges + ...
                    log(my_w) - log(1 - my_w);
            end
            
            % Accept proposal with given probability
            if log_ar > log(unifrnd(0,1))
                Theta(m, k) = theta_prop;
                Theta(k, m) = theta_prop;
                
                % Increment acceptance rate
                ar_gamma(k, m) = ar_gamma(k, m) + 1 / (burnin + nmc);
            end
            
            % Within model model
            if Theta(k, m) ~= 0
                n_within_model(k, m) = n_within_model(k, m) + 1;
                theta_prop = gamrnd(alpha_prop, beta_prop);
                Theta_prop = Theta;
                Theta_prop(k, m) = theta_prop;
                Theta_prop(m, k) = theta_prop;
                
                % Get terms that are a sum over all edges on log scale
                sum_over_edges = 0;
                for i = 1:p-1
                    for j = i+1:p
                        sum_over_edges = sum_over_edges + ...
                            log(calc_mrf_C(Theta, nu(i, j))) + 2 * ...
                            (theta_prop - Theta(m, k)) * adj(i, j, k) * adj(i, j, m) - ...
                            log(calc_mrf_C(Theta_prop, nu(i, j)));
                    end
                end
                
                % Calculate MH ratio on log scale
                log_theta_ar = (alpha - alpha_prop) * (log(theta_prop) - log(Theta(m, k))) + ...
                    (beta - beta_prop) * (Theta(m, k) - theta_prop) + sum_over_edges;
                
                % Accept proposal with given probability
                if log_theta_ar > log(unifrnd(0,1))
                    Theta(m, k) = theta_prop;
                    Theta(k, m) = theta_prop;
                    
                    % Track number of proposals accepted
                    ar_theta(k, m) = ar_theta(k, m) + 1;
                end
            end
        end
    end
    
    % Generate independent proposals for q from beta(a_prop, b_prop) density
    for i = 1:p-1
        for j = i+1:p
            q = betarnd(a_prop, b_prop);
            nu_prop = log(q) - log(1-q);
            
            % Calculate MH ratio on log scale
            % log(p(nu_prop)) - log(p(nu)) + log(q(nu)) - log(q(nu_prop))
            log_nu_ar = (nu_prop - nu(i, j)) * (sum(adj(i, j, :)) + a - a_prop) - ...
                (a + b - a_prop - b_prop) * log(1 + exp(nu_prop)) - ...
                log(calc_mrf_C(Theta, nu_prop)) + ...
                (a + b - a_prop - b_prop) * log(1 + exp(nu(i, j))) + ...
                log(calc_mrf_C(Theta, nu(i, j)));
            
            if log_nu_ar > log(unifrnd(0,1))
                nu(i, j) = nu_prop;
                nu(j, i) = nu_prop;
                ar_nu(i, j) = ar_nu(i, j) + 1 / (burnin + nmc);
            end
        end
    end
    
    % Retain values for posterior sample
    if iter > burnin
        nu_save(:, :, iter-burnin) = nu;
        Theta_save(:, :, iter-burnin) = Theta;
    end
end

% Compute acceptance rate for theta as number of proposals accepted divided
% by number of within model moves proposed
for k = 1:K-1
    for m = k+1:K
        ar_theta(k, m) = ar_theta(k, m) / n_within_model(k, m);
    end
end