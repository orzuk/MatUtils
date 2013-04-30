% Set the diagonal of a matrix to specified values
% 
function A = set_diag(A, vals)

n = length(A); % assume matrix is square 
A(1:n+1:n*n) = vals;
