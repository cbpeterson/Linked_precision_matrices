function rnorm = rMNorm(m,V,T)
% Note that this function was originally provided by Hao Wang for the
% paper Wang and Li (2012)
% rMNorm generates a qxT array  of normal draws
%  from the p-dim N(m,V)
%
p=length(m); 
rnorm = repmat(reshape(m,p,1),1,T) + chol(V)'*randn(p,T);



