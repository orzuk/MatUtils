% Find discrete and fine integrals around the unit circle
% Simulate 
% Input: 
% n - polynomial degree
% 
% Output: 
% num_bigger - # of indices with A >= R (out of 2^n)
% bigger_inds - indices of polynomials with A >= R
% C_det - determinant of corresponding cyclomatic matrix 
%
function [num_bigger bigger_inds bigger_vec C_det] = polynomial_integral(n)

R = zeros(2^n,1); A=R; C_det = zeros(2^n,1);
V = 2.*inds_to_binary(1:2^n)-1;
for k=0:(2^n-1) % loop on all possible coefficients
    k_is = k
%    v = 2.*inds_to_binary(k+1 + 2^n)-1; % set +/- 1 coefficients
%    v = v(2:end);
    v = V(k+1,:); 
    R(k+1) = 0;
    for j=0:n-1
        x = exp(2*pi*1i*j/n);
        R(k+1) = R(k+1) + abs(sum(v.*x.^(0:n-1)));
    end
    
    N = 10000; % resolution of fine integral 
    for j=0:N % Take a fine grid to approximate integral 
        x = exp(2*pi*1i*j/N);
        A(k+1) = A(k+1) + abs(sum(v.*x.^(0:n-1)));
    end    
%    A(k+1) = quadl(@(t)inner_polynomial(t, n, v), 0, 2*pi) ./ (2*pi); %  compute integral using matlab - results not stable

    C_det(k+1) = det(circulant(v)); % New: add cyclomatic matrix
    

end
R = (R ./ n); % normalize by n
A = A ./ N; % normalize by fine grid

num_bigger = sum(A >= R); 
bigger_inds = V(find(A >= R),:); 
bigger_vec = double(A >= R); 

% Internal function for computing polynomial with +/-1 coefficients 
function ret = inner_polynomial(x, n, v)

x = exp(1i.*x);
ret = zeros(size(x));
for i=0:n-1
    ret = ret + sum(v(i+1).*x.^i);
end

ret = abs(ret); %  ./ (2*pi);