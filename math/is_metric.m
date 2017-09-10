% Check if a matrix is a metric. We check that the diagonal
% is zero, that the matrix is symmetric, and that the triangle 
% inequality hold. Everything is checked within a tolerance of epsilon
%
% Input: 
% D - matrix of pairwise distances
% epsilon - tolerance
% 
% Output: 
% out_flag - flag saying if metric properties hold
% 
function out_flag = is_metric(D, epsilon) 

out_flag = 1; % assume TRUE

if(min(D(:)) < 0)  % everything must be positive
    out_flag = 0; return;
end

if(max(diag(D)) > epsilon) % check that self-distances are zero
    out_flag = 0; return;
end

if(max(max(abs(D-D'))) > epsilon) % check that matrix is symmetric
    out_flag = 0; return;
end

N = length(D);
for j=1:N
    for k=1:N
        if( min(D(:,j) + D(:,k) - D(j,k)) < -epsilon )
            out_flag = 0; return;
        end
    end
end
