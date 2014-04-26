% Compute familial risk associated with a locus with a given genetic relative risk
% 
% Input: 
% f_vec - risk allele frequency 
% grr_vec - genetic relative risk
% 
% Output: 
% lambda_s_vec - relative risk for sibling / dyzigotic twins
% lambda_mz_vec - relative risk for monozigotic twins
%
function [lambda_s_vec lambda_mz_vec ] = genetic_relative_risk_to_familial_risk(f_vec, grr_vec)

n = length(f_vec); % # of SNPs
lambda_s_vec = zeros(n,1);

for i=1:n % loop on SNPs. Each SNP represented by ONE genotypes. We do DIPLOID model 
    i_is = i; % compute lambda_s - method from Eliana
    q = f_vec(i); % get MAF
    T = [1-q+q^2/4, q-q^2/2, q^2/4; ...
        1/2-3*q/4+q^2/4, 1/2+q/2-q^2/2, q/4+q^2/4; ...
        (1-q)^2/4, 1/2-q^2/2, 1/4+q/2+q^2/4]; % transfer matrix
    D = diag( [(1-q)^2, 2*q*(1-q), q^2] ); % genotype frequencies matrix
    v = [1 grr_vec(i) grr_vec(i)*grr_vec(i)]; % assume multiplicative model
    lambda_s_vec(i) = v * D * T * v';
end
lambda_s_vec = lambda_s_vec ./ ( (1 - f_vec).^2 + ...
     2 .* f_vec .* (1-f_vec) .* grr_vec + ...
     f_vec .^ 2 .* grr_vec .^ 2).^2; % correct for prevalence (assume all are the same)

% Direct computation of lambda_s via haploid snp: 
 % lambda_s_vec2 = 0.5*((1-f_vec).*(2-f_vec)+  (2.*f_vec.*(1-f_vec)).*grr_vec + ...
%     f_vec.*(1+f_vec).*grr_vec.^2) ./ (1-f_vec+f_vec.*grr_vec).^2; 
% lambda_s_vec2=lambda_s_vec2^2
 
 lambda_mz_vec = (1-f_vec+f_vec.*grr_vec.^2) ./ (1-f_vec+f_vec.*grr_vec).^2; % simple formula - should work 
lambda_mz_vec = lambda_mz_vec.^2; % account for two alleles

%lambda_s_vec = 0.5*((1-f_vec).*(2-f_vec) +(f_vec.^2+(1-f_vec).^2).*grr_vec + ...
%    f_vec.*(1+f_vec).*grr_vec.^2) ./ (1-f_vec+f_vec.*grr_vec).^2; % simple formula - should work





