% Compute the asymptotic filtering error, when the data goes to infinity,
% and the transitions are rare
function asym_filt_err  = zeitouni_asymptotic_filtering_error(HMM_MODEL) 
M_size = length(HMM_MODEL.M);

% First determine little p: 
p_lit = sum(1-diag(M))/M_size;


% Now compute the Kullback-Leibler matrix K(i,j) :  This is not so easy for
% gaussians ! 
KL = zeros(M_size);

asym_filt_err =- p_lit * log(p_lit);