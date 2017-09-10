% Check if a matrix is 'transitive', i.e. whether it holds that: 
% A(i,j) >= 0, A(j,k) >= 0 ----> A(i,k) > 0  
%
% Input: 
% A - matrix of pairwise distances
% 
% Output: 
% tran_flag - flag saying if transitivity holds
% 
function tran_flag = is_transitive(A)

N = length(A);
tran_flag = 1; % At the beginning we assume that the matrix IS transitive
for i=1:N
    for j=1:N
        for k=1:N
            if( (A(i,j) >= 0) && (A(j,k) >= 0) && (A(i,k) < 0) )
                tran_flag = [i,j,k];
                return;
            end
        end
    end
end
            
