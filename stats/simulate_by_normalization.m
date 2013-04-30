% Simulate uniform r.v.s. normalized by their sum 
function r_vec = simulate_by_normalization(n_points)

r_vec = rand(1, n_points);

r_vec = r_vec ./ sum(r_vec); % normalize

r_vec = cumsum(r_vec); % transfer to actual values on [0,1] interval

