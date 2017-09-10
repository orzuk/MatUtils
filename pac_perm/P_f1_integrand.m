% The integrand of f1 in the saddle point
function F = P_f1_integrand( x, c_vec)
%%%%F =  (1-normcdf(x-c_vec(1))) .* normcdf(x-c_vec(2)) .* normcdf(x-c_vec(3));  % old bad ..

% New : Do the correction with density : 
F =  (1/sqrt(2*pi)) .* exp( -(x-c_vec(1)).^2./2 ) .* normcdf(x-c_vec(2)) .* normcdf(x-c_vec(3));