% Compute overlap joint distribution for linear correlations distribution 
function F = JointDensFracLinearDist(x, sigma, alpha)
ONE_SIDE = 1; TWO_SIDES = 0; 

F = (1) .* quad('lin_normcdf', -sqrt(6), sqrt(6), [], [], x, sigma, ONE_SIDE) - 1 + alpha; % Here always call with one side !!!




