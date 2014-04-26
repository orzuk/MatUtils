% Temporary calculation for joint distribution of family members liabilities

N = 1000; iters = 5000; 
X = 2.* (rand(N,iters) < 0.5) - 1; % generate first sibling
X_s = 2.* (rand(N,iters) < 0.5) - 1; % generate first sibling
X_s(1:(N/2),:) = X(1:(N/2),:); % make half genotypes equal

sigma = 0.8; % 
X_vec = sigma .* vec2column(sum(X)) ./ sqrt(N); X_s_vec = sigma .* vec2column(sum(X_s))./ sqrt(N);
E_xy = mean(X_vec .* X_s_vec)
E_x = mean(X_vec)
E_y = mean(X_s_vec)
sibling_corr = corr(X_vec, X_s_vec)

Sigma_mat = cov(X_vec, X_s_vec)
