% Compute Jacobian matrix for spherical coordinates in n dimensions
function [J, det_J] = SphericalJacobian(r, phi_vec)

n = length(phi_vec) + 1; % set dimension

J = zeros(n);

sin_vec = sin(phi_vec);
cos_vec = cos(phi_vec);

for i=1:n
    for j=1:min(i+1,n)
        J(i,j) = prod(sin_vec(1:min(i, n-1)));
        
        if(j>1) % multiply by r
            J(i,j)=r*J(i,j);
        end
        if(j == i+1) % add minus sign
            J(i,j) = -J(i,j);
        end
        
        if(i==n) % last row (special)
            if(j>1)
                J(i,j) = J(i,j) * cos_vec(j-1) / sin_vec(j-1);
            end
        else % all other rows
            if(j <= i) % don't apply for diagonal above main diagonal
                J(i,j) = J(i,j) * cos_vec(i) / sin_vec(i);
                
                if(j==1) % first column
                    
                else  % other columns
                    J(i,j) = J(i,j) * cos_vec(j-1) / sin_vec(j-1);
                end
            end
        end
        
    end
end

det_J = r^(n-1) * prod(sin_vec(1:n-2) .^ ((n-2):(-1):1));
