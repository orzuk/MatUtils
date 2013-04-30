% Perform PCA analysis specialy for genotypes
% Here we have a special 'metric' : AA = 0, AB=BA=1, BB=2
% We also add a column of '1's
function [V Y] = PachterPolytopeAnalysis(x, data_flag)

if(data_flag == 1) % special genotypes flag
    x(x == 3) = 2; % Transfer 3 to 2: 
end

x(end+1,:) = 1;  % Add one at the end for special geometry for genotypes

x = x - repmat(mean(x,2), 1, size(x,2));

% Perform S.V.D. - still doesn't work ... 
%% [U,S,V] = svd(x');
%% size(U)
%% size(S)
%% size(V)
%% Y = U'*S; % return the transformed data

% Another approach: PCA using covariance matrix
C = cov(x');
[eig_vecs eig_vals] = eig(C);
[sorted_vals sorted_inds] = sort(-diag(eig_vals));
V = eig_vecs(:,sorted_inds);
Y =  V' * x;
