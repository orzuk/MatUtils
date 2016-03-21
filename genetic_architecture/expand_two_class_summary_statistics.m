% Extract internal coding of genotype summary vec 
% 
% Input: 
% X - vector of summary statistics. Format: [Number of alleles per individual ; Number of null alleles per individual;  Number of individual carriers per allele]
%                                       or: [Number of individual carriers per allele; Number of total individual profiled per allele]
% num_individuals - number of individuals profiled
% format_flag - what is the format of X (optional) 
% 
% Output: 
% num_alleles_vec - number of alleles per individual 
% num_null_alleles_vec - number of null alleles per individual 
% num_carriers_vec - num of individuals carrying each allele
% num_total_individuals_vec - new! how many individuals were profiled for
%                             each allele (might be different for different alleles)
% L - total number of distinct alleles (includes all: null, missense, ..) 
% iters - number of different replications (simulations) 
%
function [num_alleles_vec, num_null_alleles_vec, num_carriers_vec, num_total_individuals_vec, L, iters] = ...
    expand_two_class_summary_statistics(X, num_individuals, format_flag)

iters = size(X, 2); % Input is an array of size [2*L+n] x [iters], or [2*L] x [iters]

if(~exist('format_flag', 'var') || isempty(format_flag))
    format_flag = 'individual';
end
switch format_flag
    case 'individual'
        num_alleles_vec = X(1:num_individuals,:);
        num_null_alleles_vec = X(num_individuals+1:2*num_individuals,:);
        num_carriers_vec = X(2*num_individuals+1:end,:);
        L = size(X, 1) - 2*num_individuals;
        num_total_individuals_vec = repmat(num_individuals, L, 1); % just replicate # of individuals L times
    case 'summary' % here give only cumulative number of individuals for each allele. No individual-specific information
        num_alleles_vec = []; num_null_alleles_vec = []; 
        L = size(X, 1) / 2; % each allele is repeated twice
        num_carriers_vec = X(1:L,:); % here we've got only number of people carrying each allele
        num_total_individuals_vec = X(L+1:end,:); % total number of individuals? (in case they're different!!) 
end
