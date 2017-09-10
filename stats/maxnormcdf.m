% Cumulative distribution function of maximum of N gaussians
function ret = maxnormcdf(x, N)

ret = normcdf(x).^N;
