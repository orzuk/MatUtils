% Compute average to get C_alpha from alpha 
function A = find_c_alpha(miu,sigma,prior,alpha,one_side_flag,c_a)

num_of_Gaussians=length(miu);
A=0;
for i=1:num_of_Gaussians
         A=A+prior(i)*(1-normcdf(c_a,miu(i),sigma(i))+(1-one_side_flag)*normcdf(-c_a,miu(i),sigma(i)));
end
A=A-alpha;

