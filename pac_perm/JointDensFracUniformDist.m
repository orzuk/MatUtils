% Compute overlap joint distribution for uniform correlations distribution 
function F = JointDensFracUniformDist(x, sigma, alpha)


% F = alpha - 1*alpha * normcdf( (x + 2*sqrt(3)*alpha-sqrt(3))/sigma ); % + ((sqrt(3)-x)/(2*sqrt(3))) * ...
%(normcdf( (x+2*sqrt(3)*alpha-sqrt(3))/sigma ) - normcdf( (x-sqrt(3))/sigma )) + ...
%(sigma/(2*sqrt(6*pi))) * (    exp(-(x - sqrt(3))^2/(2*sigma^2)   ) -   exp(-(x + 2*sqrt(3)*alpha - sqrt(3))^2/(2*sigma^2)   )     );
 F = (sigma/(2*sqrt(3))) * quad('normcdf', (x-sqrt(3))/sigma, (x+sqrt(3))/sigma) - 1 + alpha;
