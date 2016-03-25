% New function! plot fraction of nulls below a specific threshold
function plot_rho_finite_size(x_vec, n_samples, lambda, alpha, ...
    demographic_models_struct, i_d, new_figs_dir) % these give the population


rho_vec = frac_null_conditional_on_freq_less_f(  ...
    0.01, [], alpha, x_vec, ...
    demographic_models_struct, demographic_models_struct.model_str{i_d}) % use new function


rho_vec_finite_sample = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
    0.01, [], alpha, x_vec, ...
    demographic_models_struct, demographic_models_struct.model_str{i_d}, n_samples) % use new function - finite sample

singleton_inds = find(x_vec >= 1/(2*n_samples));

figure; semilogx(x_vec, rho_vec, 'linewidth', 2); hold on;
semilogx(x_vec(singleton_inds), rho_vec_finite_sample(singleton_inds), 'r', 'linewidth', 2);
xlabel('derived allele frequnecy, f', 'interpreter', 'latex');
ylabel('proportion of nulls, $\rho_s(f)$', 'interpreter', 'latex');
h_leg = legend({'pop. freq.', 'finite-sample'}, 1); %legend('boxoff');
set(h_leg,'xcolor',[0.8 0.8 0.8],'ycolor',[0.8 0.8 0.8]);
xlim([10^(-5) 1]); add_faint_grid(0.5);

my_saveas(gcf, fullfile(new_figs_dir, 'rho_finite_sample'), {'epsc', 'pdf', 'jpg'});
