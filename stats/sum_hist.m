% Add values in two histograms
%
% Input:
% x1 - x values for first histogram
% p1 - probablilites for first histogram
% x2 - x values for second histogram
% p2 - probablilites for second histogram
% scale_factor - resolution increase in output histogram
% normalize_flag - whether to normalize the input&output distributions (optional. Defauly: YES).
%
% Output:
% x_sum - bin values of sum histogram
% p_sum - prob. values of sum histogram
% p1 - new p1 (interpolated)
% p2 - new p2 (interpolated)
% 
function [x_sum, p_sum, p1, p2] = sum_hist(x1, p1, x2, p2, scale_factor, normalize_flag)

if(~exist('scale_factor', 'var') || isempty(scale_factor)) % set refinement level
    scale_factor = 2;
end
if(~exist('normalize_flag', 'var') || isempty(normalize_flag)) % set default normalization
    normalize_flag = 1;
end

x_sum = union(x1, x2); sum_min = min(x_sum); sum_max = max(x_sum);
new_n = round(length(x_sum) * scale_factor);
new_step = (sum_max - sum_min) / (new_n-1);
x_sum = sum_min:new_step:sum_max; % set new bins

% sum values % First generate a refinement of the distributions and bring them both to the same x vec
if(normalize_flag == 2)
    p1_integral = sum(p1); %  integral_hist(x1, p1);
    p2_integral = sum(p2); % integral_hist(x2, p2);
end
if(length(x1)>1)
    p1 = interp1(x1, p1, x_sum, 'linear', 0);
end
if(length(x2)>1)
    p2 = interp1(x2, p2, x_sum, 'linear', 0);
end

switch normalize_flag % normalize histograms to get area one (TODO: can also normalize to get area same as p1, p2 !!!!)
    case 1 % normalize to integrate to one
        p1 = normalize_hist(x_sum, p1);
        p2 = normalize_hist(x_sum, p2);
    case 2 % normalize to have same integral as original p1, p2
        p1 = p1 .* p1_integral ./ max(1/flintmax, sum(p1)); % integral_hist(x_sum, p1));
        p2 = p2 .* p2_integral ./ max(1/flintmax, sum(p2)); % integral_hist(x_sum, p2));
end % switch
p_sum = p1 + p2; % sum values in new bins

if(normalize_flag == 1) % normalize result
    p_sum = normalize_hist(x_sum, p_sum);
end
