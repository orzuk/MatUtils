% Checking new procedures for controlling the FDR 
epsilon = 0.000000001;
q = 0.05; 
num_iters = 10000; m = 2;  m0 = m; 

pvals_vec = rand(m, num_iters); % independent vars

m0_hat_amit = min(m, 2 .* sum(pvals_vec)); % the estimator of amit for m0
m0_hat_zuk = min(m, - sum(log(1-pvals_vec)));  % my estimator for m0 

% perform the FDR procedure
[sorted_pvals sort_perm] = sort(pvals_vec); % sort pvalues
tmp_mat_amit = (sorted_pvals <= repmat([1:m]' .* q, 1, num_iters) ./ repmat(m0_hat_amit, m, 1));  % perform step-up procedure
tmp_mat_zuk = (sorted_pvals <= repmat([1:m]' .* q, 1, num_iters) ./ repmat(m0_hat_zuk, m, 1));  % perform step-up procedure
tmp_mat_BH = (sorted_pvals <= repmat([1:m]' .* q, 1, num_iters) ./ m);  % perform step-up procedure
R_amit = zeros(1, num_iters); V_amit = zeros(1, num_iters);
R_amit_step_down = zeros(1, num_iters); V_amit_step_down = zeros(1, num_iters);
R_zuk = zeros(1, num_iters); V_zuk = zeros(1, num_iters);
R_zuk_step_down = zeros(1, num_iters); V_zuk_step_down = zeros(1, num_iters);
R_BH = zeros(1, num_iters); V_BH = zeros(1, num_iters);
for i=1:num_iters
    j = find(tmp_mat_amit(:,i), 1, 'last'); j_step_down = find(~tmp_mat_amit(:,i), 1, 'first');
    if(~isempty(j))
        R_amit(i) = j;
        V_amit(i) = sum(sort_perm(1:R_amit(i),i) <= m0);
    end
    if(~isempty(j_step_down))
        R_amit_step_down(i) = j_step_down-1;
        V_amit_step_down(i) = sum(sort_perm(1:R_amit_step_down(i),i) <= m0);
    end
    j = find(tmp_mat_zuk(:,i), 1, 'last'); j_step_down = find(~tmp_mat_zuk(:,i), 1, 'first');
    if(~isempty(j))
        R_zuk(i) = j;
        V_zuk(i) = sum(sort_perm(1:R_zuk(i),i) <= m0);
    end
    if(~isempty(j_step_down))
        R_zuk_step_down(i) = j_step_down-1;
        V_zuk_step_down(i) = sum(sort_perm(1:R_zuk_step_down(i),i) <= m0);
    end
    j = find(tmp_mat_BH(:,i), 1, 'last');
    if(~isempty(j))
        R_BH(i) = j;
        V_BH(i) = sum(sort_perm(1:R_BH(i),i) <= m0);
    end
end

Q_amit = V_amit ./ max(R_amit, epsilon);  Q_amit_step_down = V_amit_step_down ./ max(R_amit_step_down, epsilon); 
Q_zuk = V_zuk ./ max(R_zuk, epsilon); Q_zuk_step_down = V_zuk_step_down ./ max(R_zuk_step_down, epsilon); 
Q_BH = V_BH ./ max(R_BH, epsilon); 

FDR_amit = mean(Q_amit)
FDR_amit_step_down = mean(Q_amit_step_down)
FDR_zuk = mean(Q_zuk)
FDR_zuk_step_down = mean(Q_zuk_step_down)
FDR_BH = mean(Q_BH)

