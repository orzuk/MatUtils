% Compute mean time until absorption for allele with selection coefficient s
% The function computes the following expression:
% \int_{x_min}^{x_max} x^{weight_flag} * t_s(x) dx
% There is NO Normalization
% For equilibrium we have: t_s(f) =  (1 - e^(-S(1-f))) / [f(1-f)(1-e^(-S))]
%
% Input:
% s - selection coefficient (should be POSITIVE for deleterious alleles)
% theta - effective mutation rate (=4*N*mu)
% N - effective population size
% x_min - minimal allele frequency (default is ~0, or 1/2N)
% x_max - maximal allele frequency (default is ~1, or 1-1/2N)
% weight_flag - how do we weight different allele frequencies.
%
% Output:
% T - mean time until absorption
%
function T = absorption_time_by_selection(s, theta, N, x_min, x_max, weight_flag)


if(length(N) > 1) % allow vector input for N
    T = zeros(size(N)); 
    if(length(theta)==1)
        theta = repmat(theta, size(N)); 
    end
    for i=1:length(N)
        T(i) = absorption_time_by_selection(s, theta(i), N(i), x_min, x_max, weight_flag);
    end
    return; 
end
num_x = max(length(x_min), length(x_max));

num_s = length(s); S = 4.*N.*s; % small s to big S
if(~exist('x_min', 'var') || isempty(x_min)) % set default 
    x_min = 1/(2*N);
end
if(~exist('x_max', 'var') || isempty(x_max))
    x_max = 1-1/(2*N);
end
% two_side_flag = 0; % look at derived allele
if(~exist('weight_flag', 'var') || isempty(weight_flag))
    weight_flag = 0; % default is integrating over allele frequencies (no weights) 
end

if(x_min < x_max) % here we don't have integral, just density
    T = absorption_time_by_selection_indefinite_integral(S, theta, x_max, weight_flag) - ...
        absorption_time_by_selection_indefinite_integral(S, theta, x_min, weight_flag);
else % here we don't have integral, just density
    T = zeros(num_s, num_x);
    for i=1:num_s
        T(i,:) = theta .*  exp(allele_freq_spectrum(x_min, -S(i) ./ (4*N), N, 0, 'log', weight_flag)); %     x_hist; % here we multiply by the density
    end
end


% Internal function for computing indefinite integral !!!! 
% Input: 
% S - effective selection coefficient (S=4*N*s. should be POSITIVE for deleterious alleles)
% theta - effective mutation rate (=4*N*mu)
% x - allele frequency
% weight_flag - how do we weight different allele frequencies.
%
% Output: 
% T - Indefinite integral: int_x t_s(x) * x^{weight_flag} dx 
%
function T = absorption_time_by_selection_indefinite_integral(S, theta, x, weight_flag)

num_s = length(S); num_x = length(x); 
T = zeros(num_s,num_x);

for i=1:num_s
    switch weight_flag % which integral to compute 
        case 0 % Compute: \int t_s(x) dx
            if(S(i) ~= 0)
                T(i,:) = theta .* ( Ei(-S(i).*(x-1)) +log(x) - log(1-x) + exp(S(i)) .* ( - Ei(-S(i).*x)) ) ./ ...
                    (1 - exp(S(i))); % compute analytic solution
                if(isnan(T(i)) || isinf(T(i))) % here S is too big.
                    T(i) = -theta .* ( Ei_over_exp(-S(i).*(x-1)).*exp(x.*(-S(i))) - Ei(-S(i).*x)  );
                end
            else % neutral
                T(i,:) = theta .* log(x);
            end
        case {1, 'freq'} % average by allele frequency. Compute \int x t_s(x) dx.
            if(S(i) ~= 0)
                T(i,:) = theta .* ( Ei(-S(i).*(x-1)) - log(1-x) ) ./ (1 - exp(S(i))); % compute analytic solution
                if(isnan(T(i))) % here S is too big. Use asymptotics 
                    T(i,:) = -theta .* ( exp(-S(i).*x_max) ./ (1-x_max) ) ./ S(i);
                end
            else % integral at s=0 is equal to z
                T(i,:) = theta .* x; 
            end
        case 2 % average by allele frequency squared. Compute \int x^2 T_s(x) dx. (missing) 
            T(i,:) = ( Ei(S.*(x-1)) - log(1-x) - x + exp(S.*(x-1))./S ) ./ (1-exp(-S));
            T(i, S == 0) = x(S==0).^2./2;
            null_inds = find(isnan(T(i,:))); % here S is too big - use asymptotics.
            T(i,null_inds) = -exp(-S(null_inds) .* x(null_inds)) ./ S(null_inds).^2;

        case {-2, 'var', 'het'} % average variance explained. Compute \int x(1-x) T_s(x) dx. (DO NOT! multiply by 2). We also have an analytic solution here. (Do not multiply by theta!)
            if(S(i) > 300) % take limiting formula as -S->infinity. Why one??? should be ... 
                T(i,:) = -exp(-x.*S(i)) ./ S(i); % 0 <= f < 1 
                T(i,x==1) = -exp(-S(i)); % f=1 
            else
                if(S(i) ~= 0)
                    T(i,:) =  ( x - exp (S(i).*(1-x))./ (-S(i)) ) ./ (1-exp(S(i)));
                else % take limiting formula as S->0
                    T(i,:) = x .* (1-x./2);
                end
            end
            T(i,:) = theta .* T(i,:); % weight by theta
            
        otherwise % here take a power > 2 
            if(weight_flag > 0) % Compute indefinite integral:  \int_x x^{weight_flag} * t_s(x) dx
                if(S(i) ~= 0) 
%                    T(i,:) = 444444; % USE INTEGRATION BY PARTS !!!                    (We don't have an answer here! need a recursive formula) 
                else % S = 0
                    T(i,:) = (x.^weight_flag) ./ weight_flag; % compute integral of f^{k-1}
                end
            else % Compute integral:  \int_{x_min}^{x_max} x^{-weight_flag-1}*(1-x) * t_s(x) dx
                if(S(i) ~= 0)
                    T(i,:) = 2 .* ( (-S(i).*x).^(-weight_flag) + exp(S(i)*(1-x)) .* factorial(-weight_flag-1) .* ...
                        sum( (S(i)*f).^(0:(-weight_flag-1)) .* (-1)^(weight_flag) ./ factorial(0:(-weight_flag-1))) ) ./ ...
                        ((-S(i)).^(-weight_flag) .* (1-exp(S(i)))); % Get sum from integration by parts !!! 
                else % S = 0
                    T(i,:) = (x.^(-weight_flag) ) ./ weight_flag - ...
                        (x.^(-weight_flag-1)) ./ (weight_flag+1); % compute integral of f^{k-1}(1-f) % compute integral of f^{k-1}
                end
            end % if weight_flag > 0 (for general case) 
    end % switch weight_flag
end % loop on s



