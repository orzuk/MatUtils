% Plot comparison of selection coefficients obtained from two different sources
function plot_two_selection_comparison(HumanChimp, HumanExome, sfs_figs_dir)


% Intersect list with the list of selection coefficients from human-exome studies
[common_genes, I, J] = intersect(upper(HumanChimp.Symbol), upper(HumanExome.gene));

% Transform using quantile-normalization
[~, sort_perm] = sort(HumanExome.mis_z(J));
quant_norm_mis_z = sort(HumanChimp.s(I));
quant_norm_mis_z = quant_norm_mis_z(inv_perm(sort_perm));

figure; plot(HumanExome.mis_z(J), quant_norm_mis_z, '.')
z_vec = normcdf(-HumanExome.mis_z(J)) - 0.1 * HumanExome.mis_z(J);
z_vec = z_vec * 2*10^(-4) - 2*10^(-4);
figure; plot(HumanChimp.s(I), z_vec, '.'); hold on;
[rho, rho_pval] = corr(HumanChimp.s(I), z_vec)
xlabel('$Human-Chimp \quad \hat{s}$', 'interpreter', 'latex', 'fontsize', 14);
ylabel('$SFS \quad \hat{s}$', 'interpreter', 'latex', 'fontsize', 14);
[beta, rho2] = polyfit(HumanChimp.s(I), z_vec, 1);
pos_inds = find(HumanChimp.s(I) > (-4)*10^(-4));
[beta_pos, rho2_pos] = polyfit(HumanChimp.s(I(pos_inds)), z_vec(pos_inds), 1);
x_vec = min(HumanChimp.s(I)):0.0000001:max(HumanChimp.s(I));
plot(x_vec, x_vec*beta(1)+beta(2), 'r', 'linewidth', 2);
%plot(x_vec, x_vec*beta_pos(1)+beta_pos(2), 'r', 'linewidth', 2);
legend({'data', ['fit (\rho-' num2str(rho, 3) ')']}, 'location', 'northwest', ...
    'fontsize', 14); legend('boxoff');
my_saveas(gcf, fullfile(sfs_figs_dir, 'inter_intra_species_s_hat'), {'epsc', 'pdf', 'jpg'});


% Finally, plot the scatter of the two:
figure; hold on;
plot(HumanChimp.s(I), HumanExome.mis_z(J), '.'); % make a scatter plot of s inferred from exome data and from human-chimp comparison
xlabel('Human Chimp s estimator'); ylabel('SFS s estimator');
title(['\rho=' num2str(corr(HumanChimp.s(I), HumanExome.mis_z(J)), 3)]);


