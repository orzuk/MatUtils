% The function tries to find the 'best' curve, i.e. the curve which is as
% far as possible from the points given in the input. This curve is likely 
% to seperate two clusters. Currently we look for quadratic curves, of the 
% form : Y = a*x^2 + b*x. Output simply the values of a and b giving
% the best fit
function ab_final = AvoidPointsNewtonsMethod(ab_init, x, y, sigma, iters)

ab_current = ab_init; gamma = 0.000001;

for i=1:iters
%%    f_is = AvoidPointsGaussianFields(ab_current,x, y, sigma); 
    grad = AvoidPointsGaussianFieldsGrad(ab_current,x, y, sigma);
%%%    ab_current = ab_current - AvoidPointsGaussianFields(ab_current,x, y, sigma) .*  grad ./ (grad(1)^2+grad(2)^2);
    ab_current = ab_current + gamma .*  grad;
end

ab_final = ab_current;

