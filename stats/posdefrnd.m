% Sample a semi positive-definite matrix
% 
% Input: 
% N - matrix size
% 
% Output: 
% C - (semi) positive-definite matrix
%
function C = posdefrnd(N)

C = rand(N); C=C*C'; % get a pos-def matrix 
std_mat = sqrt(diag(1./diag(C)));
C = std_mat * C* std_mat;  % normalize to have variance one 
C=0.5*(C+C'); % just fix rounding errors
