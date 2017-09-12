% First two moments of ratio of gaussians 
function [mu sigma] = rationormstat(m1,s1,m2,s2)

mu = quadl(@(x) x.*rationormpdf(x,m1,s1,m2,s2), m1-10*(s1+s2),m1+10*(s1+s2));
sigma = sqrt(quadl(@(x) x.^2.*rationormpdf(x,m1,s1,m2,s2), mu-10*(s1+s2),mu+10*(s1+s2)) - mu^2);
