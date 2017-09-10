% Generate a sample from Ewen's distribution 
function s = ewen_sampling_formula(n, num_points,  theta)

s = chinese_restaurant_process(n, num_points, 0, theta); % A special case of chinese restuarant process 



