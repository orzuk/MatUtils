% compmute moments of |z|^p when p is standard Gaussian  
function  [mu sigma] = powernormstat(p) 

mu = 2^(p/2) * gamma((p+1)/2) / sqrt(pi); 
sigma = sqrt(2^(p) * gamma(p+1) / sqrt(pi) - mu^2);

