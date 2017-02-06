% Compute allele frequency distribution in a population
% Formula from Sawyer&Hartl paper. We use 4*N*s (not 2*N*s) since N
% is the number of individuals (not chromosomes)
%
% Input:
% x - allele frequency
% s - selection coefficient (this is little s. Negative s means deleterius and will be the standard)
% N - population size
% two_side_flag - return Derived (0) between 0 and 1 or Minor Allele Freq. (1) between 0 and 0.5
% scale_mode - return results on linear or log scale (latter to avoid underflow). Default is linear
% weight_flag - how do we weight different allele frequencies (default: no weight)
%
% Output:
% g - probability of allele frequency being below x (NOT normalized!)
%
function g = allele_freq_cumulative(x, s, N, two_side_flag, scale_mode, weight_flag)

if(~exist('two_side_flag', 'var') || isempty(two_side_flag))
    two_side_flag = 1; % default is MAF !!!
end
if(~exist('scale_mode', 'var') || isempty(scale_mode))
    scale_mode = 'linear';
end
if(~exist('weight_flag', 'var') || isempty(weight_flag))
    weight_flag = 0; % default is integral of t_s (allele-freq)
end
S = 4.*N.*s; % small s to big S
switch scale_mode
    case 'linear'
        if(two_side_flag) % collapse x and 1-x
            g = allele_freq_cumulative(x, s, N, 0, scale_mode, weight_flag) + ...
                allele_freq_cumulative(1-x, s, N, 0, scale_mode, weight_flag);
            return;
        end
        
        switch weight_flag % choose which moment to compute
            case {0} % integral of allele frequencies
                if(isscalar(S) && (S == 0))
                    g = log(x);
                else
                    g = ( log(x./(1-x)) + Ei(S.*(x-1)) - exp(-S).* ( Ei(S.*x) + pi*1i) ) ./ (1-exp(-S)); % prevent overflow for large s!
                end
                I = find(isnan(g) | isinf(g)); % get rid of nans
                if(~isempty(I)) % large S. Use asymptotics 
                    if(length(x) == 1)
                        x = repmat(x, 1, length(S));
                    end
                    g(I) = exp(S(I).*x(I)) ./ (S(I).*(x(I)-1)) + Ei(S(I).*x(I)) + pi*1i;
                end
            case 1 % multiply by f
%                g = ( Ei(S.*(x-1)) - log(1-x) ) ./  (1-exp(-S));
                g = (phi_s_integral(x, S, 1) - phi_s_integral(0, S, 1)) ./ ...
                    (phi_s_integral(1-1/(2.*N), S, 1) - phi_s_integral(0, S, 1)); % Normalized !! 
            case 2 % multiply by f^2
                g = ( Ei(S.*(x-1)) - log(1-x) - x + exp(-S.*(1-x))./S ) ./  (1-exp(-S));
            case {-2, 'var'} % multiply by f(1-f)
                g = (x - exp(-S.*(1-x))./S) ./ (1-exp(-S));
                I = find(isnan(g));
                if(~isempty(I)) % large S. Use asymptotics 
                    if(length(x) == 1)
                        x = repmat(x, 1, length(S));
                    end
                    g(I) = exp(S(I).*x(I)) ./ S(I); 
                end
        end % switch weight_flag
        
    case 'log' % compute log of allele frequency density
        if(two_side_flag) % collapse x and 1-x
            g = sum_log(allele_freq_cumulative(x, s, N, 0, scale_mode, weight_flag), ...
                allele_freq_cumulative(1-x, s, N, 0, scale_mode, weight_flag));
            return;
        end
        if(isscalar(s) && (s == 0))
            g = log(log(x));
        else
            g = log(exp(-S.*(1-x)) - 1);
            % % %             I = find(-S.*(1-x) > 300); % problematic indices (overflow)
            % % %             g(I) = -S.*(1-x(I));
            % % %             if(~var_explained_flag)
            % % %                 g = g - log(x) - log(1-x);
            % % %             end
            % % %             if(-S > 300) % assume s is a scalar
            % % %                 g = g + S;
            % % %             else
            % % %                 g = g - log(exp(-S)-1);
            % % %             end
        end % deal with general case
end % switch scale mode (linear or logarithmic)


