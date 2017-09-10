% Perform 'unique' but allow a certain errors (tolerance)
% indices should be working now ...
%
% Input:
% x - vector/matrix
% TOL - maximum tolerance for declaring 'identical'
%
% Output:
% u - unique values
% u_inds - indices: u = x(u_inds)
% v_inds - indices other-way: x = u(v_inds)
%
function [u, u_inds, v_inds] = unique_with_tol(x, TOL)

if(isvector(x))
    x = vec2row(x);
end
% First do a 'standard' unique
[x, u_inds, v_inds] = unique(x', 'rows');

[x, sort_perm] = sortrows(x); x=x';

iter=0; cont_flag = 1;
while(cont_flag == 1)
    iter = iter+1
    for i=1:size(x,2)-1
        j=i+1;
        while( (j <= size(x,2)) && (abs(x(1,i)-x(1,j)) < TOL) ) % check if we can merge
            if( max(abs(x(:,i)-x(:,j))) < TOL ) % merge i and j
                x(:,j) = x(:,i);
            end
            j=j+1;
        end
    end
    [u, u_inds2, v_inds2] = unique(x', 'rows'); u=u';
    u_inds = u_inds(u_inds2); v_inds = v_inds2(v_inds);
    if(size(u,2) == size(x,2)) % stop! no new reduction
        cont_flag = 0;
    end
    x=u;
end % while on iters

