% Test effect of directional selection and additive breeding value
% for non-additive genetic architectures
AssignGeneralConstants;
plot_flag = 0;
z_1_g_vec = -1:0.1:1; % set z ??
h_x = 0.6; % set pathway heritability
N_vec = [1:5 10]; % set number of pathways
figs_dir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/breeding/'; % save figs here

breeding_value = cell(max(N_vec),1); e_z_given_z_1_g =  cell(max(N_vec),1);
for N=N_vec % [1:3] % loop on different N's, including N=1, when the model is linear
    [breeding_value{N} e_z_given_z_1_g{N}] = ...
        compute_breeding_value_LP(h_x, N, 1, z_1_g_vec); % compute analytically (works only for linear!!!) 
end

if(plot_flag) % Plot
    figure; plot(z_1_g_vec, e_z_given_z_1_g{1});% Plot parent genotype vs. offspring phenotype for additive model
    title('parent genotype vs. offspring phenotype - breeding value, for additive model');
    xlabel('Parents genotype'); ylabel('conditional offspring phenotype shift');
    
    
    figure; hold on; % Plot parent genotype vs. offspring phenotype for non-additive LP model
    for j= 1:length(N_vec)
        N=N_vec(j);
        plot(z_1_g_vec, e_z_given_z_1_g{N}, color_vec(j));
        title(['parent genotype vs. offspring phenotype - breeding value, for LP(N' ',' ...
            num2str(h_x*100,2) '%) model']);
    end
    plot(0, 0, 'r*');
    legend([repmat('N=', length(N_vec),1) num2str(N_vec')], 4);
    xlabel('\Delta(z_{1,g})'); ylabel('\Delta(z^{offspring})');
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig_breeding_value'), 'epsc');
    
    figure; hold on;  % Plot parent genotype vs. offspring phenotype for non-additive LP model
    zero_ind = find(z_1_g_vec == 0); % index of z_1=0
    for j= 1:length(N_vec)
        N=N_vec(j);
        plot(z_1_g_vec, (e_z_given_z_1_g{N} - e_z_given_z_1_g{N}(zero_ind)) .* ...
            e_z_given_z_1_g{1}(end) ./ (e_z_given_z_1_g{N}(end)-e_z_given_z_1_g{N}(zero_ind)), color_vec(j));
        title(['parent genotype vs. offspring phenotype - breeding value, for LP(N' ',' ...
            num2str(h_x*100,2) '%) model']);
    end
    plot(0, 0, 'r*');
    legend([repmat('N=', length(N_vec),1) num2str(N_vec')], 4);
    xlabel('\Delta(z_{1,g})'); ylabel('\Delta(z^{offspring})');
    my_saveas(gcf, ...
        fullfile(figs_dir, 'supp_fig_breeding_value_normalized'), ...
        'epsc');
end % end if plot


N=1; iters = 500000; h_x = 0.0; meta_iters =3; h_shared_env =0.5*(1-h_x); % 0.25; %  0.1;
%simulate_directional_selection(iters, 2, h_x)
for i=1:meta_iters
    run_simulation_i = i
    [~, ~, ~, ~, ...
        ~, ~, ...
        ~, ~, ~, ...
        mu_s_vec(i,:) mu(i,:) mu_offspring_vec(i,:) mu_offspring(i,:)] = ...
        simulate_directional_selection(iters, N, h_x,h_shared_env, 0);  % need to insert shared environment !!! test first linear model (or 4 pathway model)
end
mean_mu = mean(mu);
mean_mu_offspring = mean(mu_offspring);
mean_mu_s_vec = mean(mu_s_vec);
mean_mu_offspring_vec = mean(mu_offspring_vec);

figs_dir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/breeding';
figure; hold on; % Plot regression of delta_mu and delta_s for different selection coefficients.
plot(mean_mu_s_vec - mean_mu, mean_mu_offspring_vec - mean_mu_offspring, '.', 'markersize', 10);
xlabel('\Delta s', 'fontsize', 14, 'fontweight', 'bold'); %  = \mu_s - \mu');
ylabel('\Delta \mu', 'fontsize', 14, 'fontweight', 'bold'); %  = \mu_{offspring} - \mu');
title([repmat(' ', 1, 80) 'supp-fig8'], 'fontsize', 16, 'fontweight', 'bold');
x_vec = min(mean_mu_s_vec-mean_mu-0.1):0.01:max(mean_mu_s_vec-mean_mu+0.1); % 0:0.1:1;
beta_offspring = polyfit(mean_mu_s_vec - mean_mu, mean_mu_offspring_vec - mean_mu_offspring, 1);
plot(x_vec, x_vec .* beta_offspring(1) + beta_offspring(2), 'r', 'linewidth', 2);
xlim([0, max(mean_mu_s_vec-mean_mu).*1.01]);
ylim([0, max(max(mean_mu_offspring_vec - mean_mu_offspring).*1.01, ...
    max(max(mean_mu_s_vec-mean_mu).*1.01 .* beta_offspring(1) + beta_offspring(2)).*1.01)])
my_saveas(gcf, fullfile(figs_dir, ...
    'response_to_selection'), {'epsc', 'pdf', 'jpg', 'fig'});



beta_selection = zeros(3,1); h_all = zeros(3,1); % vector of coefficients (?)
for i=1:20
    run_i = i
    [beta_selection(i), h_all(i), mid_parent_vec{i}, offspring_mean_vec{i} ...
        mean_selected_parents_phenotype(i) mean_selected_offspring_phenotype(i) ...
            mid_parent_qtl_vec offspring_qtl_vec selected_inds] = ...
        simulate_directional_selection(iters, 4, 0.5, 0.5*0.5, 0) % LP model P* appearing in the paper
end
mean(beta_selection); h_all = mean(h_all); 
mean_selected_parents_phenotype = mean(mean_selected_parents_phenotype);
mean_selected_offspring_phenotype = mean(mean_selected_offspring_phenotype); 

% Plot averaging
all_mid_parent_vec = cell2vec(mid_parent_vec); all_offspring_mean_vec = cell2vec(offspring_mean_vec); 
[all_mid_parent_vec sort_perm] = sort(all_mid_parent_vec); 
all_offspring_mean_vec = all_offspring_mean_vec(sort_perm);
all_offspring_mean_vec_smoothed = smooth(all_mid_parent_vec, all_offspring_mean_vec, 10);
plot_flag2=1;
if(plot_flag2)
    rand_inds = randperm(length(offspring_qtl_vec)); rand_inds = rand_inds(1:20000); 
    rand_selected_inds = intersect(selected_inds, rand_inds);
    rand_unselected_inds = setdiff(rand_inds, rand_selected_inds);
    T = 1; % selected_inds = find(mid_parent_qtl_vec_normalized > T);
    figure; hold on; x_vec = linspace(min(mid_parent_qtl_vec), max(mid_parent_qtl_vec), 1000); 
    plot(mid_parent_qtl_vec(rand_unselected_inds), offspring_qtl_vec(rand_unselected_inds), '.', 'markersize', 4); 
    plot(mid_parent_qtl_vec(rand_selected_inds), offspring_qtl_vec(rand_selected_inds), 'r.', 'markersize', 4);
    plot(all_mid_parent_vec, all_offspring_mean_vec_smoothed, 'm', 'linewidth', 2); % Plot regression of delta_mu and delta_s for different selection coefficients.
    xlabel('mid-parent phenotype'); ylabel('offspring phenotype');
    beta_fit = polyfit(all_mid_parent_vec', all_offspring_mean_vec_smoothed, 1);
    plot(x_vec, x_vec .*beta_fit(1)+beta_fit(2), 'r--', 'linewidth', 2);
    plot(x_vec, x_vec .* h_all, '--g', 'linewidth', 2); % plot line with h_all
    line([T T], [min(offspring_qtl_vec) max(offspring_qtl_vec)], 'color', 'k', 'linewidth', 3, 'linestyle', '--');
    plot(mean_selected_parents_phenotype, mean_selected_offspring_phenotype, '*k', 'linewidth', 6);
    %    text(mean_selected_parents_phenotype, mean_selected_offspring_phenotype, ...
    %        'x', 'fontsize', 20, 'fontweight', 'bold');
    legend({'unselected', 'selected', 'offspring mean', 'linear-regression', 'y=h_{all}^2 x', ...
        'selection-threshold', '\Delta \mu vs. \Delta s'}, 2);
    legend('boxoff');
    title([repmat(' ', 1, 80) 'supp-fig9'], 'fontsize', 16, 'fontweight', 'bold');
    xlim([-3,3]); ylim([min(offspring_qtl_vec(rand_inds)) max(offspring_qtl_vec(rand_inds))]); 
    my_saveas(gcf, fullfile(figs_dir, 'response_to_selection2'), {'fig', 'jpg', 'epsc'});
end


test_bivariate=0;
if(test_bivariate)% Test bi-variate distribution (what do we test here?)
    h = 0.2; iters = 10000000;
    z = randn(iters,1).*sqrt(h);
    eps = randn(iters,1).*sqrt(1-h);
    eps1 = randn(iters,1).*sqrt(1-h);
    z2 = z+eps;
    z1 = z+eps1;
    rho = corr(z, z2)
    var(z)
    var(eps)
    I = find((z2 <= 1.01) & (z2>=0.99))
    length(I)
    mean(z2(I))
    mean(z(I))
    var(z(I))
    (1-h)*h
    mean(z1(I))
    var(z1(I))
    (1-h^2)
    var(z1)
end


