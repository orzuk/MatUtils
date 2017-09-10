% Adds white gaussian noise to a vector
% 
% Input: 
% v - clean vector
% sigma - st.d. of noise
% 
% Output: 
% v_n - noisy version of vector
% 
function v_n = add_noise_to_vec(v, sigma)

z = randn(size(v,1), size(v,2)) .* sigma;
v_n = v+z;
