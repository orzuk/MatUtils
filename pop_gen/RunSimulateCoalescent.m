figs_dir = 'C:\Teaching\StatisticalEvolutionaryGenomics_52887\Docs\figs'; % replace to path for output figures

%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate for Qu 1.b. %
%%%%%%%%%%%%%%%%%%%%%%%%
n=20; theta=2; M = 5000; % set parameters (pop. size, mut. rate, #simulations)
[Mut, Trees] = SimulateCoalescent(n, theta, M);
figure; hist_density(cell2vec(Mut.ages), 500); xlabel('Age $T$',  'interpreter', 'latex', 'fontsize', 14); ylabel('Freq.', 'fontsize', 14);
title(['Empirical density of ages. $E[T] = ' num2str(Mut.age_mean, '%.3f') ', Var[T] = ' num2str(Mut.age_var, '%.3f') '$.'], 'interpreter', 'latex', 'fontsize', 14);
% Compute mean and variance
my_saveas(gcf, fullfile(figs_dir, 'coalescent_age_sim'), {'jpg'});

%%%%%%%%%%%%%%%%%%%%
% Plot for Qu 1.c. %
%%%%%%%%%%%%%%%%%%%%
t_vec = [0.1 0.5 1 2]; ages_vec = cell2vec(Mut.ages);  counts_vec = cell2vec(Mut.counts); figure;
for i=1:length(t_vec)
    good_inds = find( abs( ages_vec - t_vec(i) ) < 0.1 ); % find mutations with this age
    allele_freq = hist(counts_vec(good_inds), 1:(n-1));
    subplot(2,2,i); hold on; bar(1:(n-1), allele_freq ./ length(good_inds));
    plot( 1:(n-1), 1 ./ (sum(1./(1:(n-1))) .* (1:(n-1))), 'r', 'linewidth', 2);
    if(i==2)
        legend({'Empirical', '$\frac{1}{i H_{n-1}}$'}, 'interpreter', 'latex', 'fontsize', 14); legend('boxoff');
    end
    ylabel('Freq.');
    if(i >= 3)
        xlabel('Derived Allele Count');
    end
    title(['$t= ' num2str(t_vec(i)) '\pm 0.1$'], 'interpreter', 'latex', 'fontsize', 14);
end
suptitle('Empirical distribution of derived allele counts for different allelic ages');
my_saveas(gcf, fullfile(figs_dir, 'coalescent_count_given_age_sim'), {'jpg'});

%%%%%%%%%%%%%%%%%%%%
% Plot for Qu 1.d. %
%%%%%%%%%%%%%%%%%%%%
mean_age_vec = zeros(n-1,1); var_age_vec = zeros(n-1,1);
for k=1:(n-1)
    mean_age_vec(k) = mean(ages_vec(counts_vec == k));
    var_age_vec(k) = var(ages_vec(counts_vec == k));
end
figure; errorbar( 1:(n-1), mean_age_vec, sqrt(var_age_vec), 'linewidth', 2);
xlabel('Derived Allele Count', 'fontsize', 14); ylabel('Age', 'fontsize', 14);
title('$E[T] \pm st.d.(T)$ as function of derived allele count', 'interpreter', 'latex', 'fontsize', 14);
my_saveas(gcf, fullfile(figs_dir, 'coalescent_age_given_count_sim'), {'jpg'});
