% Compute the mean and variance of z given two loci
% The input h_x_l is ASSUMING ONE LIABILITY !!!!
%
% Input:
% N - total # of liabilities
% K - # liabilities needed to exceed threshold
% mu - prevalence
% h_x - total narrow-sense heritability of the trait assuming LT model 
% h_x_l - heritability explained by x_1 for the trait's liability assuming LT model
%
% Output:
% z_mu - prevalence
% z_var - (average) variance of z given x_1,x_2. Integrate over all (x_1,x_2) values
% z_var_one - ??? 
% 
function [z_mu z_var z_var_one] = z_stats_given_two_x_MLT(N, K, mu, h_x, h_x_l)


if(~exist('h_x_l', 'var') || isempty(h_x_l))
    h_x_l = heritability_scale_change_MLT(N*h_x, K, N, mu, 'MLT'); % Get the heritabiity of each liability in MLT model
end
options = optimset('tolx', 0.0000000000000000000001); % increase optimization tolerance
mu_l = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu)), 0, 1, options); % find mu_l that keeps the prevalence
x_mu_l = norminv(1-mu_l); % set threshold

if(N == 1) %  A single liability
    z_mu = quadl(@(x)normpdf(x).*(1-normcdf( (x_mu_l - x.*sqrt(2.*h_x_l)) ./ sqrt(1-2*h_x_l))),-10,10);
    z_var = quadl(@(x)normpdf(x).* ...
        (1-normcdf( (x_mu_l - x.*sqrt(2.*h_x_l)) ./ sqrt(1-2*h_x_l))).* ...
        normcdf( (x_mu_l - x.*sqrt(2.*h_x_l)) ./ sqrt(1-2*h_x_l)),-10,10);
else % at least two liabilities. Loci are in different liabilities
    z_mu = 1 - (1-mu_l).^(N-2) .* ...
        quadl(@(x)normpdf(x).*normcdf((x_mu_l-x.*sqrt(h_x_l))./sqrt(1-h_x_l)),-10,10) .* ...
        quadl(@(x)normpdf(x).*normcdf((x_mu_l-x.*sqrt(h_x_l))./sqrt(1-h_x_l)),-10,10);
    z_var = 1-z_mu - (1-mu_l).^(2*(N-2)) .* ...
        quadl(@(x)normpdf(x).*normcdf((x_mu_l-x.*sqrt(h_x_l))./sqrt(1-h_x_l)).^2,-10,10) .^2; %...
    %    quadl(@(x)normpdf(x).*normcdf((x_mu_l-x.*sqrt(h_x_l))./sqrt(1-h_x_l)).^2,-10,10);
    z_var_one = 1-z_mu - (1-mu_l).^(2*(N-1)) .* ...
        quadl(@(x)normpdf(x).*normcdf((x_mu_l-x.*sqrt(h_x_l))./sqrt(1-h_x_l)).^2,-10,10); %...
    %     z_var_one2 =  (1-mu_l).^(N-1) .* quadl(@(x)normpdf(x).* ...
    %         (1-(1-mu_l).^(N-1)*normcdf( (x_mu_l - x.*sqrt(h_x_l)) ./ sqrt(1-h_x_l))).* ...
    %         normcdf( (x_mu_l - x.*sqrt(h_x_l)) ./ sqrt(1-h_x_l)),-10,10);
end

