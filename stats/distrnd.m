% Simulate data from any distribution (like normrnd but for general dists.)
% 
% Input: 
% dist_x - x vector of distribution 
% dist_p - probabilities at each x
% m - width of output 
% n - length of output
%
% Output: 
% data - a data matrix with values taken from dist
%
function data = distrnd(dist_x, dist_p, m, n, varargin)
if(~exist('m', 'var'))
    m=1;
end
if(~exist('n', 'var'))
    n=1;
end
k = length(dist_x); % number of bins in the distribution 
data = rand(m*n,1); [data sort_perm] = sort(data);  % take data from a uniform distribution and sort it
dist_p = cumsum(dist_p ./ sum(dist_p)); % get a cumulative sum 
cur_ind = 1; 
for i=1:k
   new_ind = find(data(cur_ind:end) > dist_p(i), 1); 
   if(~isempty(new_ind))
        data(cur_ind:new_ind+cur_ind-2) = dist_x(i);
        cur_ind=cur_ind+new_ind-1;
   else % this means that we've passed the whole distribution - finish the loop 
       data(cur_ind:end) = dist_x(i); 
       break;
   end
end
data = reshape(data(sort_perm), m, n); % sort back to order and reshape 





