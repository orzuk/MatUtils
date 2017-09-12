% Multiply two elements in GF(2^m).
%
% Input:
% m - power of 2
% x - first element
% y - second element
% p - irreducible polynomial
%
% Ouptut:
% r - result: x*y mod p
%
function r = mult_gf2m(m, x, y, p)


r = 0;
for i=1:m
    r = bitxor(r, bitget(x, i)*2^(i-1)*y);
    for j=2*m:-1:m+1
        if(bitget(r, j))
            r = bitxor(r, p*2^(j-m-1));
        end
    end
end



