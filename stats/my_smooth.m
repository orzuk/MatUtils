% Performs moving sum (default)/average
% It is needed since it appears as if matlab smoothing
% has a bug for very large vectors of singles for some reason.
%
% Input:
% v - a vector of numbers
% k - a smoothing window size
% ave_flag - flag saying if to compute moving sum (0, default) or average(1)
% shift - distance between two consecutive windows
% do_edges - flag saying if to compute average also at edges (1) or not (0, default)
%
% Output:
% v_smoothed - a vector of length n-k+1 containing moving average of v
%
function v_smoothed = my_smooth(v, k, ave_flag, shift, do_edges, varargin)

if(~exist('k', 'var') || isempty(k)) % set default window size
    k = 5;
end
k = min(k, length(v)); % can't have smooth window longer than vec
if(~exist('ave_flag', 'var') || isempty(ave_flag))
    ave_flag = 0;
end
if(~exist('shift', 'var') || isempty(shift))
    shift = 1;
end
if(~exist('do_edges', 'var') || isempty(do_edges))
    do_edges = 0;
end
row_flag = is_row(v); % check if its a row vector
v = vec2column(v); % make sure it's a column vector)
n = length(v); k = double(k);


if(do_edges)
    k = k-1+mod(k,2); % force k to be odd - problem: this doesn't allow to smooth with even numbers!!!     
    m = floor(n / shift); % size of output
else
    m = floor((n-k) / shift) + 1; % size of output
end
if(isa(v, 'single'))
    v_smoothed = zeros(m,1, 'single'); % do smoothing 'manually'
else
    v_smoothed = zeros(m,1); % do smoothing 'manually'
end
if(do_edges)
    v_smoothed = k .* smooth(double(v), k); % simply call the matlab smooth function (doesn't work for large singles)
    v_smoothed(1:floor(k/2)) = v_smoothed(ceil(k/2));
    v_smoothed(end-floor(k/2):end) = v_smoothed(end-ceil(k/2)); % set edges 
    if(isa(v, 'single'))
        v_smoothed = single(v_smoothed);
    end
else
    for i=1:k
        v_smoothed = v_smoothed + v(i:shift:i+shift*(m-1));
    end
end
if(row_flag)
    v_smoothed = vec2row(v_smoothed);
end
if(ave_flag)
    v_smoothed = v_smoothed ./ k;
end
