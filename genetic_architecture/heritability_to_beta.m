% Compute effect size from heritability for a QTL
%
% Input:
% h - heritability value of QTL
% f_vec - minor allele frequencies
%
% Output:
% beta - effect size (regression)
%
function beta = heritability_to_beta(h, f_vec, SNP_TYPE)

if(~exist('SNP_TYPE', 'var') || isempty(SNP_TYPE))
    SNP_TYPE = 'binary';
end
switch SNP_TYPE
    case 'diploid' % multiply by two
        beta = sqrt(h ./ (2.*f_vec.*(1-f_vec)));
    case 'binary' 
        beta = sqrt(h ./ (f_vec.*(1-f_vec)));
end

