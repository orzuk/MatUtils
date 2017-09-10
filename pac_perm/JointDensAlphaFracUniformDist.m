% Compute overlap joint alpha distribution for uniform correlations distribution 
function F = JointDensAlphaFracUniformDist(y, sigma, alpha)

%%%sigma = 1;
F = 1 - normcdf( (y + sqrt(3))/sigma ) + ((sqrt(3)-y)/(2*sqrt(3))) * ...
(normcdf( (y+sqrt(3))/sigma ) - normcdf( (y-sqrt(3))/sigma )) + ...
(sigma/(2*sqrt(6*pi))) * (    exp(-(y - sqrt(3))^2/(2*sigma^2)   ) -   exp(-(y + sqrt(3))^2/(2*sigma^2)   )     ) - alpha;

