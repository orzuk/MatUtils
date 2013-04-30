% The integrand of f0 in the saddle point
function F = P_f0_integrand( x, c_vec)
%%%F =  normcdf(x-c_vec(1)) .* (     normcdf(x-c_vec(2)) .*(1-normcdf(x-c_vec(3))) + (1-normcdf(x-c_vec(2))) .* normcdf(x-c_vec(3))); %% old bad

% New : Do the correction with density : 
F =  (1/sqrt(2*pi)) .* normcdf(x-c_vec(1)) .* ...
    (     normcdf(x-c_vec(2)) .* exp( -(x-c_vec(3)).^2/2 ) + exp( -(x-c_vec(2)).^2/2 ) .* normcdf(x-c_vec(3))  );
