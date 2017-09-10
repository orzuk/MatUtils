% Randomize an mXn matrix such that each column sums to 1
function A = rand_normalized(m, n)

A = rand(m,n); 
A = A ./ repmat(sum(A), m, 1); 



