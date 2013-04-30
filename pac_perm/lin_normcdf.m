% Compute multiplication of linear and normal cdf 
function pp = lin_normcdf(c, z_alpha, sigma, one_side_flag)
if(one_side_flag)
    pp = (1/sqrt(6) - abs(c)/6).*normcdf((z_alpha-c)./sigma);
else
     pp = 2*(1/sqrt(6) - abs(c)/6).*(normcdf((z_alpha-c)./sigma) - normcdf((-z_alpha-c)./sigma));
end