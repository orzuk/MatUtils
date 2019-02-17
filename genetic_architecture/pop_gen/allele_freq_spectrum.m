% Compute allele frequency distribution in a population using analytic approximation
% Formula from Sawyer&Hartl paper. We use S=4*N*s (not 2*N*s) since N
% is the number of individuals (not chromosomes)
% The formula is: g(f) = t_s(f) =  [1 - e^(-S(1-f))] / [f(1-f)(1-e^(-S))].
% This is density of the time being spent at frequency, such that expected
% time spent between frequencies [a,b] is \int_a^b g(f) df
% Note: This is HALF (1/2) of the corret value !! (a factor 2 is missing here!).
% This function gives approximately same answer as allele_freq_spectrum_numeric if dividing by a factor of 2N/2 = N
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
% g - g(i) is mean-time-spent density of allele frequency being x(i) (NOT normalized!). 
%
function g = allele_freq_spectrum(x, s, N, two_side_flag, scale_mode, var_explained_flag)

if(~exist('two_side_flag', 'var') || isempty(two_side_flag))
    two_side_flag = 1; % default is Minor Allele Frequency !!!
end
if(~exist('scale_mode', 'var') || isempty(scale_mode))
    scale_mode = 'linear';
end
if(~exist('var_explained_flag', 'var') || isempty(var_explained_flag))
    var_explained_flag = 0; % default is allele-freq. (not variance explained) !!!
end
S = 4.*N.*s; % small s to big S
switch scale_mode
    case 'linear'
        % NEW! just take log part
        g = exp(allele_freq_spectrum(x, s, N, two_side_flag, 'log', var_explained_flag)); 
    case 'log' % compute log of allele frequency density
        if(two_side_flag) % collapse x and 1-x
            g = sum_log(allele_freq_spectrum(x, s, N, 0, scale_mode, var_explained_flag), ...
                allele_freq_spectrum(1-x, s, N, 0, scale_mode, var_explained_flag));
            return;
        end
        if(isscalar(s) && (s == 0))
            if(~var_explained_flag)
                g = -log(x);
            else
                g = log(1-x);
            end
        else
            g = log(exp(-S.*(1-x)) - 1);
            I = find(-S.*(1-x) > 300); % problematic indices (overflow)
            g(I) = -S.*(1-x(I));
            if(~var_explained_flag)
                g = g - log(x) - log(1-x);
            end
            if(isscalar(S) && (-S > 300)) % assume s is a scalar
                g = g + S;
            else
                g = g - log(exp(-S)-1);
            end
        end % deal with general case
end


