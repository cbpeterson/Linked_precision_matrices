clear;
%%%%%
addpath './MCMC_code';
addpath './Helper_functions';
addpath './Case_study/Code';

% Y1 = importHC("log_HC.csv");
cur_filename = strcat('./Case_study/Input/fake_hpHC.csv');
Y1 = readtable(cur_filename);
Y1 = table2array(Y1);
%%%%%%%
S1 = Y1' * Y1;

cur_filename = strcat('./Case_study/Input/fake_HC.csv');
Y2 = readtable(cur_filename);
Y2 = table2array(Y2);
%%%%%%%
S2 = Y2' * Y2;

    
cur_filename = strcat('./Case_study/Input/fake_MCI.csv');
Y3 = readtable(cur_filename);
Y3 = table2array(Y3);
%%%%%%%
S3 = Y3' * Y3;
                                                                                                                                                                                                                                                    
cur_filename = strcat('./Case_study/Input/fake_AD.csv');
Y4 = readtable(cur_filename);
Y4 = table2array(Y4);
S4 = Y4' * Y4;


data = cat(3, S1, S2, S3 ,S4);
N = [size(Y1,1),size(Y2,1),size(Y3,1),size(Y4,1)];

PPI = zeros(100,100,4); MEDIAN = zeros(100,100,4);

for si = 1:4
% Y = data(:,:,si);
% S = Y' * Y;
S = data(:,:,si);
n = N(si);
p = 100;
    
lambda = 1; pii = 2/(p-1);
h = 100; v0 = 0.1; v1 = 15;
V0 = v0*ones(p); V1 = v1*ones(p); burnin = 10000; nmc = 10000;
[Sig_save,C_save,Z_save] = BayesGGM_SSVS_FixedV0V1(S,n,S/n,V0,V1,lambda,pii,burnin,nmc);
    

PPI(:,:,si) = mean(Z_save,3);
MEDIAN(:,:,si) = mean(Z_save,3)>.5;

end

%% write files to be used in later analysis %%
csvwrite('./Case_study/Output/SingleGraph_hpHC.csv',PPI(:,:,1))
csvwrite('./Case_study/Output/SingleGraph_HC.csv',PPI(:,:,2))
csvwrite('./Case_study/Output/SingleGraph_MCI.csv',PPI(:,:,3))
csvwrite('./Case_study/Output/SingleGraph_AD.csv',PPI(:,:,4))

