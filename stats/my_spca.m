% Crude sparse-pca. Compute PCA and threshold smallest loadings
% Input: 
% X - data matrix
% s - number of non-zero loadings for each PC 
% Output: 
% W - the sparse principal components 
% 
function W = my_spca(X, s)

[~,~,W] = svd(X); % Compute PCA directions
for i=1:size(W,1)
   [~, P] = sort(abs(W(:,i))); 
   W(P(1:(size(W,2)-s)), i) = 0;
   W(:,i) = W(:,i) ./ sqrt( sum(W(:,i).*W(:,i)) ); % normalize
end
