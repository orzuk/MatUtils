% Test the function hypergeometric_for_many_sets


N = 2000; % universe size
t_vec = [140 311 800 55 90 200 110 500 660 167 800 123]; % Many sets sizes - good gaussian approximation
% t_vec = [140 311 800 ]; % Few sets. Not so great gaussian approximation
K = length(t_vec); % Number of sets

bin_mat = zeros(K,N); % the matrix of zeros and ones representing occurances
iters = 10000;
num_bins = 100; 


m_1 = zeros(K,iters); m_2 = m_1;
s_1 = cell(K,K);  s_2 = s_1;

s = zeros(iters,1); over_pval = s; under_pval = s; z_score = s;
for i=1:iters
    for j=1:K      % Sample occurances matrix
        bin_mat(j,:) = rand_nchoosek(N, t_vec(j));
    end
    m_1(:,i) = bin_mat(:,11); m_2(:,i) = bin_mat(:,2); % copy first two rows
    sets_intersection_mat = bin_mat * bin_mat';
    s(i) = ( sum(sum(sets_intersection_mat)) - sum(diag(sets_intersection_mat)) ) / 2; % get total intersection size
    [over_pval(i) under_pval(i) z_score(i) mu_s sigma_s] = hypergeometric_for_many_sets(N, t_vec, sets_intersection_mat);
end

% Compute 4 types of variance/covariance empirically
for i=1:K
    for j=1:K
        s_1{i,j} = m_1(i,:) .* m_1(j,:);
        s_2{i,j} = m_2(i,:) .* m_2(j,:);
    end
end
emp_sigma = zeros(4,1);

for i=1:K
    for j=1:K
        if(j ~= i)
            emp_sigma(1) = emp_sigma(1) + var(s_1{i,j});
            emp_sigma(2) = emp_sigma(2) + my_cov(s_1{i,j}',s_2{i,j}');
            for jj = 1:K % 3 and 4
                if( (jj ~= i) && (jj ~= j))
                    emp_sigma(3) = emp_sigma(3) + my_cov(s_1{i,j}',s_1{i,jj}');
                    emp_sigma(4) = emp_sigma(4) + my_cov(s_1{i,j}',s_2{i,jj}');
                end
            end
        end
    end
end

% Multiply by combinatorial constants
emp_sigma(1) = (N/2) * emp_sigma(1); % divide by 2 since we took i \neq j and not i<j 
emp_sigma(2) = (N*(N-1)/4) * emp_sigma(2); % divide by 2 since we took i \neq j and not i<j
emp_sigma(3) = N * emp_sigma(3); 
emp_sigma(4) = (N*(N-1)/2) * emp_sigma(4); 


figure; hist(over_pval, num_bins); title('pvals dist.');  xlabel('p-val'); ylabel('freq.');
figure; hold on; hist(z_score, num_bins); title('Z scores dist.');  xlabel('Z score'); ylabel('freq.');

figure; hold on; h = hist_density(s, num_bins); title('sum-of-intersects dist.');  xlabel('s'); ylabel('freq.');
mu_samp = mean(s); sigma_samp = std(s);
line([mu_s mu_s], [0 max(h)], 'color', 'r'); % computed mean line
line([mu_samp mu_samp], [0 max(h)], 'color', 'r', 'LineStyle', ':'); % empirical mean line
line([mu_s+sigma_s mu_s+sigma_s], [0 max(h)], 'color', 'g'); % computed +std line
line([mu_s-sigma_s mu_s-sigma_s], [0 max(h)], 'color', 'g'); % computed -std line
line([mu_samp+sigma_samp mu_samp+sigma_samp], [0 max(h)], 'color', 'g', 'LineStyle', ':'); % empirical +std line
line([mu_samp-sigma_samp mu_samp-sigma_samp], [0 max(h)], 'color', 'g', 'LineStyle', ':'); % empirical -std line
legend('hist', 'comp. mean', 'emp. mean', 'comp. std.', '', 'emp. std.', '');

% Old debugging stuff
% sigma_is = sigma'
% emp_sigma_is = emp_sigma'
% emp_sigma_s = sqrt(sum(emp_sigma))
% sigma_samp_is = sigma_samp
% mean_diff_why_pos = mu_s - mu_samp
