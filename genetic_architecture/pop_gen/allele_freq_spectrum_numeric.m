% Compute allele frequency distribution in a population using Markov Chain numeric calculation
% Note: This is Normalized to one !! 
%
% Input:
% x - allele frequency
% s - selection coefficient (this is little s). Should be NEGATIVE for delterious alleles (so fitness is 1+s)
% N - population size
% two_side_flag - return Derived (0) between 0 and 1 or Minor Allele Freq. (1) between 0 and 0.5 (Default: 1)
% scale_mode - return results on linear or log scale (latter to avoid underflow). Default is linear
% var_explained_flag - if this is 'ON', we compute var explained distribution instead
%                      (multiply by x*(1-x)). Default is 'OFF'.
%
% Output:
% g - probability of allele frequency being x (NOT normalized!)
%
function g = allele_freq_spectrum_numeric(x, s, N, two_side_flag, scale_mode, var_explained_flag)

if(~exist('two_side_flag', 'var') || isempty(two_side_flag))
    two_side_flag = 1; % default is Minor Allele Frequency !!!
end
if(~exist('scale_mode', 'var') || isempty(scale_mode))
    scale_mode = 'linear';
end
if(~exist('var_explained_flag', 'var') || isempty(var_explained_flag))
    var_explained_flag = 0; % default is allele-freq. (not variance explained) !!!
end

switch scale_mode
    case 'log' 
        g = log(allele_freq_spectrum_numeric(x, s, N, two_side_flag, 'linear', var_explained_flag));        
    case 'linear'
        if(two_side_flag) % collapse x and 1-x
            g = allele_freq_spectrum_numeric(x, s, N, 0, scale_mode, var_explained_flag) + ...
                allele_freq_spectrum_numeric(1-x, s, N, 0, scale_mode, var_explained_flag);
            return;
        end
        M = FisherWright_ComputeMarkovMatrix(N, s, 'exact', 1); % compute Markov matrix 
        mu=1;         M(1,2)=mu; M(1,1)=1-mu; M(end,2)=mu; M(end,end)=1-mu; % add small mutation to make process ergodic 
        g = markov_chain_stationary_dist(M); g = g(2:(end-1)) ./ sum( g(2:(end-1)) ); % get stationary distribution 
end % switch scale 


