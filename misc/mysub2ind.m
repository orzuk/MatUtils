
% Compute sub2ind in one vector. Inverse of myind2sub 
%
% Input: 
% L - vector of lengths
% d - number of dimensions
% ix_vec - vector of indices for each dimension 
% 
% Output: 
% out - index 
% 
function out = mysub2ind(L, d, ix_vec)

if(isscalar(L))
    L = repmat(L, [1 d]); % dimension array for a d-dimension array L long on each side
end
c = num2cell(ix_vec); % cell([1 d]);  % dynamically sized varargout
out = sub2ind(L, c{:}); 


