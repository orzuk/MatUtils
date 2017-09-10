% Draw effect size uniformly at random U[0,1] 
% and plot the resulting multiplicative lambda_s
function simulate_rand_grr(iters, output_file)

%n_vec = [250:250:1000];
n_vec = 10000; % [250:250:1000];
lambda_vec = 1 + 1 ./ sqrt(n_vec);
mu = 0.01; % disease prevalence
for i=1:length(n_vec) % loop on different values of N's
    
    f_vec = rand(n_vec(i), iters); % randomize the minor allele frequencies
    %    f_vec(:) = 0.1; % set fixed effect
    lambda_s_mult = zeros(iters,1);
    tmp_lambda_vec = mat_into_vec(repmat(lambda_vec(i), n_vec(i),2));
    for j=1:iters
        do_iter = j
        tmp_f_vec = mat_into_vec(repmat(f_vec(:,j), 1, 2)');
        %        lambda_s_mult(j) = prod(tmp_lambda_vec);
        [lambda_s_vec lambda_s_add lambda_mz_add h_add V_add lambda_s_mult(j)] = ...
            genetic_relative_risk_to_heritability(tmp_f_vec, tmp_lambda_vec, mu);
    end
    lambda_s_mult_mean(i) = mean(lambda_s_mult);
    lambda_s_mult_std(i) = std(lambda_s_mult);
end


figure; hold on; errorbar(n_vec, lambda_s_mult_mean, lambda_s_mult_std, '.');
plot(n_vec, repmat(exp(1/12)*2-1, length(n_vec), 1), 'k-'); % take square for diploid snps
title('multiplicative \lambda_s for uniformly distributed allele-frequencies with effect size 1+1/\sqrt{n}');
xlabel('N'); ylabel('\lambda_s');
my_saveas(gcf, output_file, {'jpg', 'epsc', 'pdf'});  % save figure
