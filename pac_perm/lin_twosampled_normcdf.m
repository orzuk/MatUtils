% Compute multiplication of linear and normal cdf for two sampled distributions
function pp = lin_twosampled_normcdf(c, z_alpha, sigma, one_side_flag)
if(one_side_flag)
    pp = (1/sqrt(6) - abs(c)/6).* (1-normcdf((z_alpha-c)./sigma)).^2;
else
     pp = 2*(1/sqrt(6) - abs(c)/6).*(1-normcdf((z_alpha-c)./sigma)).*(1-normcdf((z_alpha-c)./sigma) + normcdf((-z_alpha-c)./sigma));
end