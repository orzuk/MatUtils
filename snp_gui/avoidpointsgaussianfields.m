% This function gives the value of the gaussian fields when we have
% parabole parameters a and b and data x,y
function C = AvoidPointsGaussianFields(ab,x, y, sigma) 
a = ab(1); b = ab(2);

C = -sum( exp( -(y-a*x.^2-b.*x).^2./(2*sigma^2) ) );







