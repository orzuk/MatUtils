% Compute q*q multipliation table for GF(2^m)
% 
% p - the primitive polynomial 
% 
function T = compute_mult_tab(m, p)

q = 2^m;
T = zeros(q);

for i=1:q
    for j=1:q
        T(i,j) = mult_gf2m(m, i-1, j-1, p);
    end
end
