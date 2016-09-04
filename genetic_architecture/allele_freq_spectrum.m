% Compute allele frequency distribution in a population
% Formula from Sawyer&Hartl paper. We use S=4*N*s (not 2*N*s) since N
% is the number of individuals (not chromosomes)
% The formula is: g(f) = t_s(f) =  [1 - e^(-S(1-f))] / [f(1-f)(1-e^(-S))]
%
% Input:
% x - allele frequency
% s - selection coefficient (this is little s). Should be NEGATIVE for delterious alleles (so fitness is 1+s)
% N - population size
% two_side_flag - return Derived (0) between 0 and 1 or Minor Allele Freq. (1) between 0 and 0.5
% scale_mode - return results on linear or log scale (latter to avoid underflow). Default is linear
% var_explained_flag - if this is 'ON', we compute var explained distribution instead
%                      (multiply by x*(1-x)). Default is 'OFF'.
%
% Output:
% g - probability of allele frequency being x (NOT normalized!)
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
        if(two_side_flag) % collapse x and 1-x
            g = allele_freq_spectrum(x, s, N, 0, scale_mode, var_explained_flag) + ...
                allele_freq_spectrum(1-x, s, N, 0, scale_mode, var_explained_flag);
            return;
        end
        %        g = (exp(-S.*(1-x)) - 1) ./ (x.*(1-x).*(exp(-S)-1)); % problem: overflow for large s!
                
        if(isscalar(S) && (S == 0))
            g = 1 ./ x;
        else
            g = (exp(S.*x) - exp(S)) ./ (x.*(1-x).*(1-exp(S))); % prevent overflow for large s!
        end        
        
        %   if(two_side_flag)
        %       g = g + (exp(-4.*N.*s.*(x)) - 1) ./ (x.*(1-x).*(exp(-4.*N.*s)-1));
        %   end
        if(var_explained_flag)
            g = g .* x .* (1-x);
        end
        
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


