% This shows the joint probability distribution. It depends on the original
% correlation function. 
function F = MixJointDensFrac(y, sig, x_alpha, one_side_flag,true_corr_flag,prior,miu,dist_std)

num_of_Gaussians=length(prior);
Y=0;
for i=1:num_of_Gaussians
    Y=Y+prior(i)*(1./(dist_std(i)*sqrt(2.*pi))) .* exp(-(y-miu(i)).*(y-miu(i))./(2*dist_std(i)^2));
end
if(true_corr_flag)
     F = Y .* (1 - normcdf( (x_alpha-y)./sig )+(1-one_side_flag)*normcdf( (-x_alpha-y)./sig ));
else
      F = Y .* (1 - normcdf( (x_alpha-y)./sig )+(1-one_side_flag)*normcdf( (-x_alpha-y)./sig )).^2;
end 
