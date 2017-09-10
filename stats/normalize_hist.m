% Normalize a histogram to have sum one 
% 
% Input: 
% x - x values
% p - counts/probabilities
% standardize_flag - make distribution have mean zero and st.d. one (default is false)
% 
% Output: 
% p_norm - normalized probabilities
% x_norm - normalized x vector (only when standardizing)
% mu - mean of standardized distribution
% sigma - std of standardized distribution
% normalization_factor - by how much did we divide to normalize (need to add)
% 
function [p_norm, x_norm, mu, sigma] = normalize_hist(x, p, standardize_flag)

if(max(p) == 0) % avoid zeros (assume only positive probs.)
    p(:) = 1/flintmax; 
end
if(length(p) > 1)
    p_norm = p ./ integral_hist(x, p); % max(integral_hist(x, p), 1/bitmax); 
else
    p_norm = 1; 
end

if(~exist('standardize_flag', 'var') || isempty(standardize_flag))
    standardize_flag = 0; 
end
mu = mean_hist(x, p_norm); sigma = std_hist(x, p_norm);
if(standardize_flag)
    x_norm = (x - mu) ./ sigma; % normalize bin sizes
else
    x_norm = x;
end
