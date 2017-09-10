% Compute mean and st.d. of maximum of N standard Gaussians
function [mu sigma] = maxnormstat(N)

mu = quadl(@(x)maxnormpdf(x,N).*x, -10, 10);
sigma = sqrt(quadl(@(x)maxnormpdf(x,N).*x.^2, -10, 10) - mu.^2);

