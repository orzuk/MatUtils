% Take the bitand of all elements in the vector
function A = bitand_vec(X)

if(~isvector(X)) % enable a matrix too
    A = zeros(size(X,1),1);
    for i=1:size(X,1)
        A(i) = bitand_vec(X(i,:));
    end
else
    X = unique(X);
    A = X(1);
    for i=2:length(X)
        A = bitand(A, X(i));
    end
end