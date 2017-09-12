% Normalize a two-dimensional histogram to have sum one
%
% Input:
% x - x values
% y - y values
% p - counts/probabilities matrix
% standardize_flag - make distribution have mean zero and st.d. one (default is false)
%
% Output:
% p_norm - normalized probabilities matrix
% x_norm - normalized x vector (only when standardizing)
% y_norm - normalized y vector (only when standardizing)
% mu - mean of standardized distribution
% sigma - std of standardized distribution
% normalization_factor - by how much did we divide to normalize (need to add)
%
function [p_norm x_norm y_norm mu sigma] = normalize_hist2d(x, y, p, standardize_flag)

if(max(p) == 0) % avoid zeros (assume only positive probs.)
    p(:) = 1/bitmax;
end
if(length(p) > 1)
    p_norm = p ./ integral_hist2d(x, y, p); % max(integral_hist(x, p), 1/bitmax);
else
    p_norm = 1;
end

% Part below not working yet
% if(~exist('standardize_flag', 'var') || isempty(standardize_flag))
%     standardize_flag = 0;
% end
% mu = mean_hist2d(x, y, p_norm); sigma = std_hist2d(x, y, p_norm);
% if(standardize_flag)
%     x_norm = x - mu(1); x_norm = x_norm ./ sigma(1); % normalize bin sizes
%     y_norm = y - mu(2); y_norm = y_norm ./ sigma(2); % normalize bin sizes   
% else
%     x_norm = x; y_norm = y; 
% end
% 
% 
