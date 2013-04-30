% Computes all types of a given length 
% 
% Input: 
% k - ??
% m - length
% count_couples - ??
% 
% Output: 
% types_vec - a vector of all types
% 
function types_vec = TypesEnumerate(k, m, count_couples)

% first generate all possible vectors of length k from alphabet of size k
temp_vecs = zeros(k, m^k);

for i=1:k
    base_vector = reshape(repmat([0:m-1], m^(i-1), 1), 1, m^i);
    temp_vecs(i,:) =  repmat(base_vector, 1, m^(k-i));
end

types_vec = zeros(m, m^k);
for i=1:m
    types_vec(i,:) = sum(temp_vecs == (i-1));
end

if(count_couples) % here ALL types are probably different ... 
    types_mat = zeros(k-1,m,m,m^k); % good only for low k and m
    for i=1:k % first position 
        for j=i+1:k % second position 
            for s1=1:m % first symbol
                for s2=1:m  % second symbol
                    f = find((temp_vecs(i,:) == s1-1) & (temp_vecs(j,:) == s2-1));
                    types_mat(j-i,s1,s2,f) = types_mat(j-i,s1,s2,f)+1;
                end
            end
        end
    end
    types_vec = reshape(types_mat, (k-1)*m^2, m^k);
end
types_vec = unique(types_vec', 'rows');

