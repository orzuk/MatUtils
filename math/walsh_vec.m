% Generate a walsh/hadamard vector of size N 
function w = walsh_vec(N)

w = [0 1]; 
for i=2:N
    w = [w 1-w];
end

       
