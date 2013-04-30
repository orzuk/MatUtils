% Multiply two polynomials
% 
% The input: 
% coeff1 - coefficient vector of first polynomial
% coeff2 - coefficient vector of second polynomial
% 
% The output: 
% coeff_vec - coefficient vector of output polynomial
%
function coeff_vec = mult_pol(coeff1, coeff2)

coeff1_sparse = sparse(coeff1);
coeff2_sparse = sparse(coeff2);

[i1,j1,s1] = find(coeff1_sparse);
[i2,j2,s2] = find(coeff2_sparse);
T = max(max(j1)-1, max(j2)-1); % max degree
num_values1 = max(size(j1)); % % % num_values2 = max(size(j2));   % old version 
coeff_vec = zeros(1, 2*T+1); % 0,...,2*T

% Now do multiplication (heaviest part ...)
for i = 1:num_values1 
    coeff_vec(j1(i)-1 + j2) = coeff_vec(j1(i)-1 + j2) + s1(i)*s2;  
end

% old version 
% % % % Now do multiplication (heaviest part ...)
% % % for i = 1:num_values1
% % %     for j = 1:num_values2
% % %         coeff_vec(j1(i)-1 + j2(j)-1 +1) = coeff_vec(j1(i)-1 + j2(j)-1 +1) + s1(i)*s2(j);
% % %     end
% % % end
