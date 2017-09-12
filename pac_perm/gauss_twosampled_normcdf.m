% Compute overlap for two samples from gaussian distribution 
function pp = gauss_twosampled_normcdf(c, z_alpha, sigma, one_side_flag)
% if(one_side_flag)
%     pp = (1/sqrt(2*pi)).* exp(-c.*c./2).*(1-normcdf((z_alpha-c)./sigma)).^2;
% else
%      pp = (2/sqrt(2*pi)).* exp(-c.*c./2).*(1-normcdf((z_alpha-c)./sigma)+ normcdf((-z_alpha-c)./sigma)).^2;  %(1-normcdf((z_alpha-c)./sigma)).*
% end


if(one_side_flag)
    pp = (1/sqrt(2*pi)).* exp(-c.*c./2).*(1-normcdf((z_alpha-c)./sigma)).^2;
else
        
     pp = (1/sqrt(2*pi)).* exp(-c.*c./2).*((1-normcdf((z_alpha-c)./sigma) + normcdf((-z_alpha-c)./sigma)).^2);   % Take the sqr
%%%%%%     pp = (2/sqrt(2*pi)).* exp(-c.*c./2).*(1-normcdf((z_alpha-c)./sigma)).*(1-normcdf((z_alpha-c)./sigma) + normcdf((-z_alpha-c)./sigma));  %(1-normcdf((z_alpha-c)./sigma)).*
end