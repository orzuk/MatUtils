% Test how not to divide the set of number
function g = test_factorial(N)

g = N;
for k=1:N-1
    x = nchoosek(N,k); g = gcd(x,g);     
end

