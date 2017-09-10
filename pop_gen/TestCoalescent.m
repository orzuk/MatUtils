% Test coalescent
n=3; theta=10; M = 100000;
[MutAnal, TreesAnal] = AnalyticCoalescent(n, theta);

[Mut, Trees] = SimulateCoalescent(n, theta, M); 
Mut.levels_freq


mean_a = 0;
for i=1:n % use scalar version to allow optimziation for vector theta
    mean_a = mean_a + theta ./ (theta + i-1);
end
mean_alleles_analytic = mean_a 
mean_alleles = mean(Mut.n_alleles)
mean_sites = mean(Mut.num)



good_inds = find(length_cell(Mut.counts) > 0); 
for i=1:M
    Mut.freqs{i} = hist(Mut.counts{i}, 1:(n-1)); Mut.freqs{i} = Mut.freqs{i}' ./ sum(Mut.freqs{i}); 
end
CCC = cell2num(Mut.freqs);
mean(CCC(:,good_inds), 2)

Theta_hat = zeros(M, 3); % get matrix of esitmators
for i=vec2row(good_inds) % 1:M
    if(mod(i, 1000) == 0)
        get_G = i
    end    
    G_mat = genotype_mat_from_indices(Trees.Topology(i,:), Mut.levels{i}, Mut.branches{i}); % generate genotype matrix
    Theta_hat(i,:) = EstimateThetaCoalescent(G_mat);
end

mean(Theta_hat) % estimate theta 
MSE = mean((Theta_hat - theta).^ 2)


figs_dir = 'C:\Users\oorzu\Google Drive\HUJI\Teaching\StatisticalEvolutionaryGenomics_52887\Docs\figs';

%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate for Qu 1.b. %  
%%%%%%%%%%%%%%%%%%%%%%%%
n=20; theta=2; M = 50000; % set parameters 

[Mut2, Trees2] = SimulateCoalescent(n, theta, M); 
[MutAnal, TreesAnal] = AnalyticCoalescent(n, theta);

figure; hist_density(cell2vec(Mut.ages), 500); xlabel('Age $T$',  'interpreter', 'latex', 'fontsize', 14); ylabel('Freq.', 'fontsize', 14); 
title(['Empirical density of ages. $E[T] = ' num2str(Mut.age_mean, '%.3f') ', Var[T] = ' num2str(Mut.age_var, '%.3f') '$.'], 'interpreter', 'latex', 'fontsize', 14); 
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
        xlabel('Dericed Allele Count'); 
    end
    title(['$t= ' num2str(t_vec(i)) '\pm 0.1$'], 'interpreter', 'latex', 'fontsize', 14); 
    suptitle('Empirical distribution of derived allele counts for different allelic ages'); 
end
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

