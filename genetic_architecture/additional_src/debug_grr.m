% function debug_grr()

iters = 1000000;

mu=0.1; x_mu = norminv(1-mu); 
h=0.05;
f=0.2; 

x = randn(iters,1) .* sqrt(1-h);

g1 = randn(iters,1) .* sqrt(h);  %here geneotype is Gaussian
g2 = ((rand(iters,1) < f)-f) * sqrt(h / (f*(1-f))); 

mean(g1)
mean(g2)

z1 = (x + g1) > x_mu;
z2 = (x + g2) > x_mu; 

gg2 = g2 /sqrt(h / (f*(1-f))) +f;
grr2 = (sum(z2 .* gg2) / sum(gg2)) / (sum(z2 .*(1-gg2)) / sum(1-gg2))

[grrrr grr_analytic] = heritability_to_genetic_relative_risk(h, 'liability', f, mu)
