% AssignStatsConstants()

ROC = 0; FDR = 1;  % How to display results
APARAM_INTERSECT = 2; APARAM_DIFF = 3; MIX_GAUSS = 4; UPPER_BOUND_EPS = 5; % methods for separating curves 
FDR_INV = 6; APARAM_CUM_DIFF = 7; 

diff_flag_vec = cell(10,1); % vector of diff_flags 
diff_flag_vec{FDR} = 'FDR';
diff_flag_vec{APARAM_INTERSECT} = 'APARAM_INTERSECT';
diff_flag_vec{APARAM_DIFF} = 'APARAM_DIFF';
diff_flag_vec{MIX_GAUSS} = 'MIX_GAUSS';
diff_flag_vec{UPPER_BOUND_EPS} = 'UPPER_BOUND_EPS';
diff_flag_vec{FDR} = 'FDR';
diff_flag_vec{FDR_INV} = 'FDR_INV';
diff_flag_vec{APARAM_CUM_DIFF} = 'APARAM_CUM_DIFF';

% different similarity/dissimilarity metrics
EUCLIDIAN = 0;
LOGLIKE = 1;
PEARSON = 2;
DOTPROD = 3;
HAMMING = 4; % number of identical elements
EDIT_DIST = 5; % This the hamming edit distance (how much do insertion/deletion/mismatach cost?
KULLBACK_LEIBLER = 6; % relative entropy distance 
JENSEN_SHANNON = 7; % symmetrized Kullback-Leibler
COVARIANCE = 8; % covariance of two r.v.s

LODS = 0; LOGODDS = 0; LOGRATIO = 0; MINUS_LOGRATIO = -1; % log-ratio (all the same thing) 
BAYESIAN = 0; ML = 1; PARSIMONY = 2; % phylogenetic methods 
AIC = 0; BIC = 1; MDL = 2; % Information criteria 

GAUSS_SMOOTH_SIGMA = 12; % default smoothing parameter for Gaussian smoothing 

MAD_CONST = 1/norminv(0.75); % converts from robust mad estimator to standard sigma estimator
MAD_CONST_SQR = MAD_CONST*MAD_CONST;

BINARY=0; DISCRETE=1; CONTINUOUS = 2; PROBS = 3; % variables types
