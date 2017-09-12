% Compare two sampled correlation with uniform distribution 
function pp = uni_twosampled_normcdf(c, z_alpha, sigma, one_side_flag)
if(one_side_flag)
    pp = (1/(2*sqrt(3))).* (1-normcdf((z_alpha-c)./sigma)).^2;
else
     pp = (1/sqrt(3)).*(1-normcdf((z_alpha-c)./sigma)).*(1-normcdf((z_alpha-c)./sigma) + normcdf((-z_alpha-c)./sigma));
end