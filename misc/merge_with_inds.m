% Merge two lists and keep the indices in the original list 
function [z inds] = merge_with_inds(x, y, top_k, varargin)

[z inds] = sort([x y]); 

if(exist('top_k', 'var')) % take only maximal 
    z = z(end-top_k+1:end); inds = inds(end-top_k+1:end); 
end









