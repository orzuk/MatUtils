%  A function for computing the mean of a large single array that
%  workarounds a matlab bug, by simply sampling fewer values  
function m = mean_single(x, iters, varargin)

% Just sample many numbers and take their mean
samp_inds = unique(unidrnd(length(x), 1000000,1)); 

m = mean(x(samp_inds)); 


% if(~exist('iters', 'var'))
%     iters = 4;
% end
% m = mean(x); 
% for i=1:iters
%    m = mean(x-m) + m;  
% end
