% Compute the mean of the variance explained curve for a given s 
% 
% Input: 
% s - selective coefficient
% N - effective population size
% var_flag - how to compute mean (count alleles, count varienace explained or consider allele frequencies)
% x_min - minimal value of x's possible 
%
% Output: 
% mu_x - mean value of the allele frequency x, provided 
% 
function mu_x = mean_x_given_s(s, N, var_flag, x_min)

S=4*N.*s;

if(~exist('var_flag', 'var') || isempty(var_flag)) % default is to compute mean of variance explained
    var_flag = 1; 
end
if(~exist('x_min', 'var') || isempty(x_min))
    x_min = 1 ./ (2.*N);
end

% if(-S > 300)
%     mu_x = 1./S.^2; % The result of integrating x*phi_s(x) * x*(1-x) dx
% else
% end


switch var_flag
    case 1 % mean of var. explained distribution 
%        mu_x = ((S-1).^2+1-2.*exp(-S)) ./ (2.*S.^2 .* (1-exp(-S))); % The result of integrating x*phi_s(x) * x*(1-x) dx
        
        mu_x = ((S-1).^2+1-2.*exp(-S)) ./ (2.*S .* (S-1+exp(-S)));
        I = find(-S > 300); % Asymptotics for very large (absolute) S
        mu_x(I) = -1 ./ S(I);
        J = find(-S < 0.0001); % Asymptotics for very small (absolute) S
        mu_x(J) = 1/3 + S(J) ./ 36 - S(J).^2 ./ 540;  
    case 0 % mean of allele frequency distribution 
        EulerGamma = -psi(1); % Euler's constant
        mu_x = log(1-x_min) - Ei(S.*(x_min-1)) + EulerGamma + log(-S); % numerator
        mu_x = mu_x ./ (mu_x - log(x_min) + exp(-S).*(Ei(S.*x_min)-Ei(S))); % denominator
end

