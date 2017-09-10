% calc the marginal of the probability P on the variables in I
function [P_C] = collapse_prob(P,I)

n = log2(size(P,2)); % Get the dimension
n_c = length(I);

if(n_c < n)
    I_bin = sum(2 .^ (I-1));
    I_inv_bin = sum(2 .^ (setdiff(1:n,I)-1));
%     P_C = zeros(size(P,1),2^n_c);
    IndVec = stretch_to_indexes(reshape( repmat(0:2^n_c-1,2^(n-n_c),1), 2^n,1), I_bin) + ...
        stretch_to_indexes(repmat(0:2^(n-n_c)-1,1,2^n_c)', I_inv_bin);
    P_C = reshape(sum(reshape(P(:,IndVec+1)',  2^(n-n_c), size(P,1)*2^n_c)),  2^n_c, size(P,1))';
else
    P_C = P;
end