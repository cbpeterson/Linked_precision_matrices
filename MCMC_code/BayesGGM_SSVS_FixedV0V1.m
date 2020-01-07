function [Sig_save,C_save,Z_save] = BayesGGM_SSVS_FixedV0V1(S,n,Sig,V0,V1,lambda,pii,burnin,nmc)

%% BayesGGM_SSVS_FixedV0V1 fits SSSL (Wang 2014) to a concentration graph models
%% Input:
%%   S: p x p cross product matrix:  Y*Y' if Y is a p x n matrix of random sample
%%   n: sample size
%%   Sig: initial guess of Sigma
%%   V0:  p x p matrix of the  small variance  components;  V0(i,j) is the small variance of the precision element C(i,j)
%%   V1:  p x p matrix of the  large variance  components;  V0(i,j) is the small variance of the precision element C(i,j)
%%   lambda: hyperparameter for the diagonal element Usually set to be 1;
%%   pii: prior marginal edge "inclusion probability"
%%   burnin: # of discarded samples
%%   nmc: # of saved samples

%% Output:
%%  Sig_save: p x p x nmc  saved covariance matrices
%%  C_save: p x p x nmc saved concentration matrices
%%  Z_save: p x p x nmc saved graph indicator matrices

C = inv(Sig);
p = size(S,1); 


 


C_save = zeros(p,p,nmc); Sig_save = C_save; Z_save = C_save; 
Z =  ones(p);



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


pii_RB = zeros(p);
pii_mat = zeros(p);


%grid1 = [0:0.05:100];
%grid2 = [-60:0.05:60];

%%logd_offdiag_save = zeros(nmc,length(grid2));
%%logd_diag_save  = zeros(nmc,length(grid1));

for iter = 1: burnin+nmc    
            
       if(mod(iter,3000)==0)
        fprintf('iter = %d nedge = %d \n',iter, (sum(Z(:))-p)/2);
       end
    
%%% sample Sig and C = inv(Sig)        
    for i = 1:p


      ind_noi = ind_noi_all(:,i);
 
       
      tau_temp = tau(ind_noi,i);
       
      Sig11 = Sig(ind_noi,ind_noi); Sig12 = Sig(ind_noi,i);
      
      invC11 = Sig11 - Sig12*Sig12'/Sig(i,i);
      
      Ci = (S(i,i)+lambda)*invC11+diag(1./tau_temp);
  
      Ci = (Ci+Ci')./2;       
      Ci_chol = chol(Ci);    
      mu_i = -Ci_chol\(Ci_chol'\S(ind_noi,i));
        beta = mu_i+ Ci_chol\randn(p-1,1);
        
        
        C(ind_noi,i) = beta;
        C(i,ind_noi) = beta;
        
        a_gam = 0.5*n+1;
        b_gam = (S(i,i)+lambda)*0.5;
        gam = gamrnd(a_gam,1/b_gam);
        
        c = beta'*invC11*beta;
        C(i,i) = gam+c;
    

        
        
        
        %% Below updating Covariance matrix according to one-column change of precision matrix
        invC11beta = invC11*beta;
        
        Sig(ind_noi,ind_noi) = invC11+invC11beta*invC11beta'/gam;
        Sig12 = -invC11beta/gam;
        Sig(ind_noi,i) = Sig12;
        Sig(i,ind_noi) = Sig12';
        Sig(i,i) = 1/gam;

        
        
              
         v0 = V0(ind_noi,i);
         v1 = V1(ind_noi,i);
                
        w1 = -0.5*log(v0) -0.5*beta.^2./v0+log(1-pii);
        w2 = -0.5*log(v1) -0.5*beta.^2./v1+log(pii);
        
        w_max = max([w1,w2],[],2);

        w = exp(w2-w_max)./sum(exp([w1,w2]-repmat(w_max,1,2)),2);

        z = (rand(p-1,1)<w);
        
        
        v = v0;
        v(z) = v1(z);
        
        
        
      
        
        pii_mat(ind_noi,i) = w;
        
        tau(ind_noi,i) = v;        
        tau(i,ind_noi) = v;

        Z(ind_noi,i) = z;
        Z(i,ind_noi) = z;

        
        
    end

    
    

      
       if iter >burnin   
            pii_RB = pii_RB + pii_mat/nmc; 
            Sig_save(:,:,iter-burnin) = Sig; 
            C_save(:,:,iter-burnin) = C;
            Z_save(:,:,iter-burnin) = Z;            

            
        %    logd_offdiag_save(iter-burnin,:) = logd_offdiag;
        %    logd_diag_save(iter-burnin,:) = logd_diag;

       end



end

     %       logd_max = max(logd_offdiag_save);
     %        d_offdiag = exp(logd_max).*mean(exp(logd_offdiag_save-repmat(logd_max,nmc,1)));


     %       logd_max = max(logd_diag_save);
     %       d_diag = exp(logd_max).*mean(exp(logd_diag_save-repmat(logd_max,nmc,1)));

        
