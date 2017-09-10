% Here mu_l, x_mu_l refer to a SINGLE liability. h_x refers to heritability
% on a liability (out of many in MLT)
% Input:
% x1 - value of first liability
% x2 - value of second liability
% N - number of liabilities 
% h_x - heritability in each liability
% mu_l - 'prevalence' of each liability (not used) 
% x_mu_l - threshold of each liability
% 
% Output: 
% z_mu - expected z given these two liabilities 
% 
function z_mu = z_expected_given_two_x_MLT(x1, x2, N, h_x, mu_l, x_mu_l)

if(N == 1)
    z_mu = 1 - normcdf((x_mu_l - sqrt(h_x).*(x1+x2))./sqrt(1-2.*h_x));
else
    z_mu = 1 - (1-mu_l).^(N-2) .* ...
        normcdf((x_mu_l - x1.*sqrt(h_x))./sqrt(1-h_x)) .* ...
        normcdf((x_mu_l - x2.*sqrt(h_x))./sqrt(1-h_x));
end


