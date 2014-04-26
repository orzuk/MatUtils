% Compute mean time until absorbtion for allele with selection coefficient s
% The function computes the following expression:
% \int_{x_min}^{x_max} x^{weight_flag} * t_s(x) dx
% There is NO Normalization
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
% T - mean time until absorbtion
%
function T = absorbtion_time_by_selection(s, theta, N, x_min, x_max, weight_flag)

num_s = length(s); S = 4.*N.*s; % small s to big S
T = zeros(num_s,1);
if(~exist('x_min', 'var') || isempty(x_min))
    x_min = 1/(2*N);
end
if(~exist('x_max', 'var') || isempty(x_max))
    x_max = 1-1/(2*N);
end
x_vec =  x_min:(1/(2*N)):x_max; % (1-1/(2*N)); %   (1:(2*N-1)) ./ (2*N);
two_side_flag = 0; % look at derived allele
if(~exist('weight_flag', 'var') || isempty(weight_flag))
    weight_flag = 0; % default is integrating over allele frequencies
end
for i=1:num_s
    switch weight_flag
        case 0 % Compute: \int T_s(x) dx
            %            x_hist = exp(allele_freq_spectrum(x_vec, -s(i), N, two_side_flag, 'log')); % We can replace this by an analytic expression!!
            if(S(i) ~= 0)
                T(i) = theta .* ( Ei(-S(i)*(x_max-1)) +log(x_max) - log(1-x_max) - Ei(-S(i)*(x_min-1)) -log(x_min) + log(1-x_min) + ...
                    exp(S(i)) .* (Ei(-S(i)*x_min) - Ei(-S(i)*x_max)) ) ./ ...
                    (1 - exp(S(i))); % compute analytic solution
                if(isnan(T(i))) % here S is too big.
                    T(i) = theta .* -( Ei_over_exp(-S(i)*(x_max-1)) - Ei_over_exp(-S(i)*(x_min-1)) + ...
                        Ei(-S(i)*x_min) - Ei(-S(i)*x_max)  );
                    %                     exp(-S(i)*x_min) ./ (x_min.*(1-x_min)) - exp(-S(i)*x_max) ./ (x_max.*(1-x_max)) ) ./ S(i);
                end
            else % neutral
                T(i) = theta * (log(x_max)-log(x_min));
            end
        case {1, 'freq'} % average by allele frequency. Compute \int x t_s(x) dx.
            %            x_hist = exp(allele_freq_spectrum(x_vec, -s(i), N, two_side_flag, 'log'));
            %            x_hist = x_hist .* x_vec;
            if(S(i) ~= 0)
                T(i) = theta .* ( Ei(-S(i)*(x_max-1)) - log(1-x_max) - Ei(-S(i)*(x_min-1)) + log(1-x_min) ) ./ ...
                    (1 - exp(S(i))); % compute analytic solution
                if(isnan(T(i))) % here S is too big. Use asymptotics 
                    T(i) = theta .* (exp(-S(i)*x_min) ./ (1-x_min) - exp(-S(i)*x_max) ./ (1-x_max) ) ./ S(i);
                end
            else % integral at s=0 is equal to z
                T(i) = theta * (x_max-x_min);
            end
        case 2 % average by allele frequency squared. Compute \int x^2 T_s(x) dx.
        case {-2, 'var'} % average variance explained. Compute \int x(1-x) T_s(x) dx. We also have an analytic solution here
            %            x_hist = exp(allele_freq_spectrum(x_vec, -s(i), N, two_side_flag, 'log', 1));
            if(S(i) ~= 0)
                T(i) =  ( (x_max - exp (S(i).*(1-x_max))./ (-S(i)) ) - (x_min - exp (S(i).*(1-x_min))./ (-S(i)) ) ) ./ ...
                    (1-exp(S(i)));
            else % take limiting formula as S->0
                T(i) = x_max .* (1-x_max./2) - x_min .* (1-x_min./2);
            end
    end
    if(length(x_vec) > 1) % perform integral
        %         switch weight_flag
        %             case {-2, 'var'} % use histogram. Why also on '1'?
        %                 T(i) = theta * integral_hist(x_vec, x_hist);
        %         end
    else
        if(x_min >= x_max)
            T(i) = theta * x_hist; % here we multiply by the density
        end
    end
    %    T2(i) = theta * integral_hist(x_vec, x_hist);
end % loop on s
