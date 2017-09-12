% Take the bitxor of all elements in the vector
function A = bitxor_vec(X)
X = unique(X); 
A = X(1);
for i=2:length(X)
    A = bitxor(A, X(i));
end
