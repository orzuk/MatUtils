% Compute mean of the the offspring's phenotype given the genotype of the parent
% 
% Input: 
% h_x - the heritability in each pathway
% N - number of pathways
% k - number of pathways needed to exceed the threshold 
% z_1_g 
%
% Output: 
% breeding_value - ??
% e_z_given_z_1_g - ?? 
%
function [breeding_value e_z_given_z_1_g] = compute_breeding_value_LP(h_x, N, k, z_1_g)

TOL = 0.00000000000001;
[mu_z sigma_z] = maxnormstat(N);
plot_density = 0; 

num_z = length(z_1_g);
for i=1:num_z
    compute_z = i
    e_z_given_z_1_g(i) = quadl(@(t) t.*maxnormcdf(t, N-2) .* ...
        ( normcdf(t) .* normpdf((t-h_x.*z_1_g(i)./2)./(1-h_x/2)) ./ (1-h_x/2) + ...
        (N-1).*normpdf(t).*normcdf((t-h_x.*z_1_g(i)./2)./(1-h_x/2)) ), ... % integrand
        mu_z-5.*sigma_z, mu_z+5.*sigma_z, TOL);
    if(z_1_g(i)==0) % plot central density
        t_vec = -15:0.01:15;
        t_norm_vec = (t_vec-h_x.*z_1_g(i)./2)./(1-h_x/2);
        y_vec = maxnormcdf(t_vec, N-2) .* ...
            ( normcdf(t_vec) .* normpdf(t_norm_vec) ./ (1-h_x/2) + ...
            (N-1).*normpdf(t_vec).*normcdf(t_norm_vec) );
        y_cumulative_vec = maxnormcdf(t_vec, N-1) .* ...
            normcdf(t_norm_vec);
        if(plot_density)
            figure;
            plot(t_vec, y_vec); title(['N=' num2str(N)]);
            figure;
            plot(t_vec, y_cumulative_vec); title(['N=' num2str(N)]);
        end
    end
end
e_z_given_z_1_g = (e_z_given_z_1_g-mu_z) ./ sigma_z; % normalize trait
breeding_value = 2.*e_z_given_z_1_g ./ z_1_g;

