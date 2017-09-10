% Compute Bartlett type correction for paper by Porteous 1985
% n - # deg. freedom
% p - size 
% 
% Output: 
% est_mean - estimator of E [ log |Q_MLE| ]
%
function est_mean = Porteous_cov_correction(n, p)

k_vec = n+1 - (1:p); 

est_mean = sum(log(k_vec) - 1./k_vec - 1./(3.*k_vec.^2));

