% Build a binary matrix representing the reed-solomon work
%
% Input:
% q - field order (=2^m)
% deg - degree of polynomial
% p_irreducible - irreducible polynomial representing the field
%
% Output:
% layered_code_mat - a matrix representing the code: q layers, q rows per layer, and q^(deg+1) columns
%
function layered_code_mat = build_layered_reed_solomon_code(q, deg, p_irreducible)

m = log2(q); % power of 2 (field is GF[2^m])
layered_code_mat = zeros(q^2, q^(deg+1));

for i = 1:q % loop on layer
    run_layer = i
    for j=1:q % loop on element in layer
        run_row = i*q+j
        for k = 1:q^(deg+1) % loop over polynomials
            coeff_vec = my_dec2base(k-1, q, deg+1); % find the individual coefficients for each polynomial
            p_x = 0; x_pow = 1;
            for jj = 0:deg
                p_x = bitxor(p_x, mult_gf2m(m, coeff_vec(jj+1), x_pow, p_irreducible));
                x_pow = mult_gf2m(m, x_pow, i-1, p_irreducible);
            end
            if(p_x == j-1) % Evaluate polynomial on x(i) and compare to y(j)
                layered_code_mat( (i-1)*q+j, k) = 1;
            end
        end
    end
end
