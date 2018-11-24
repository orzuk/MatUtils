% Unique making Nans into one value 
function [U, I, J] = unique_nan(x, flag)
if ~exist('flag','var')
    [U, I, J] = unique(x);
else
    [U, I, J] = unique(x, flag);
end
I(isnan(U(1:end-1))) = [];
U(isnan(U(1:end-1))) = [];
J(isnan(x)) = length(U); % assume last 

