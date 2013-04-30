% We compute the mean of X for various situations

time_elapsed = cputime;
N=4; p=0.5;
num_edges = N*(N-1)/2;
x_mean = zeros(1, num_edges+1); x_std = zeros(1, num_edges+1);
i=1;


[x_mean x_std p_x_positive_exact p_x_positive_2nd_moment perms_corrs] = calc_x_aligned_moments(N, p, [0:num_edges]./num_edges);

% % % % for alpha=([0:num_edges]./num_edges)
% % % %     alpha
% % % %     [x_mean(i) x_std(i)] = calc_x_aligned_moments(N, p, alpha);
% % % %     i=i+1;
% % % % end


figure; hold on; errorbar([0:num_edges], x_mean./factorial(N), x_std./sqrt(factorial(N))); 
plot([0:num_edges], p_x_positive_2nd_moment, 'r'); plot([0:num_edges], p_x_positive_exact, 'g');
title(['X mean for various \alpha vals, N=' num2str(N) ' p=' num2str(p)]); xlabel('\alpha'); ylabel('x mean');
legend('X moments', 'P(X>0) lowerbound', 'P(X>0)');


figure; hold on; plot(perms_corrs(end-1,:), perms_corrs(3,:), '.'); xlabel('Perms hamming'); ylabel('X correlation'); 
title('Relatios between permutation distance and correlation of alignment'); 

figure; hold on; plot(perms_corrs(end,:), perms_corrs(3,:), '.'); xlabel('Perms swaps'); ylabel('X correlation'); 
title('Relatios between permutation distance and correlation of alignment'); 

time_elapsed = cputime - time_elapsed


% Now plot for each 