% Density of maximum of N standard Gaussians
function ret = maxnormpdf(x, N)

ret = N .* normpdf(x) .* normcdf(x).^(N-1);

