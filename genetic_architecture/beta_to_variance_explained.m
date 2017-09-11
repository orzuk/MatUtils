% Convert the regression coefficient beta into units of variance explained
%
% Input:
% beta - coefficient of regression: Y = beta*X + alpha
% f - allele frequency
% V - variance of the trait
% SNP_TYPE - do we look at diploid SNP (real) or binary (easy modelling)
% n_cases - sample size (optional) used for correction 
%
% Output:
% V_explained - fraction of trait's variance explained by current genotype
%
function V_explained = beta_to_variance_explained(beta, f, V, SNP_TYPE, n_cases)

if(~exist('SNP_TYPE', 'var') || isempty(SNP_TYPE))
    SNP_TYPE = 'binary';
end
if(~exist('V', 'var') || isempty(V)) % default variance noramlzied to one
    V = 1;
end
switch SNP_TYPE
    case 'binary'
        V_explained = f.*(1-f) .* beta.^2 ./ V; % normlize both genotype and phenotype to have st.d.=1. Assume just 0,1
    case 'diploid'
        V_explained = 2.*f.*(1-f) .* beta.^2 ./ V; % here assume HW equilibrium and 0,1,2
end

if(exist('n_cases', 'var')) % New: apply finite sample correction to get an unbiased estimator 
    V_explained = V_explained - 1 ./ n_cases; % subtract 1/n
end
