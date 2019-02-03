% Find closest demography to a given demography
function [I_min, demographic_dist] = get_closest_demography(D, D_target)

demographic_dist = zeros(D.num_params, 1);
N_vec = demographic_parameters_to_n_vec(D_target);
for i=1:D.num_params  % loop on all models of first
    demographic_dist(i) = demographic_models_distance(N_vec, demographic_parameters_to_n_vec(D, i)); % find which one is most similar to the TRUE model
end
[~, I_min] = min(demographic_dist);
    

