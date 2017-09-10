% Plot comparison of selection coefficients obtained from two different sources
% 
% Input: 
% HumanChimp - structure representing Human-Chinpanzee differences
% HumanExome - structure representing allele frequencies in humans 
% sfs_figs_dir - where to save output figures 
% 
function plot_two_selection_comparison(HumanChimp, HumanExome, sfs_figs_dir)


% Intersect list with the list of selection coefficients from human-exome studies
[common_genes, I, J] = intersect(upper(HumanChimp.Symbol), upper(HumanExome.gene));

% here get rid of problematic genes 
if(isfield(HumanChimp, 'good_inds'))
    [I, I_G] = intersect(I, HumanChimp.good_inds);
    J = J(I_G); 
end

% Transform using quantile-normalization
if(isfield(HumanExome, 'mis_z'))
    quant_norm = 1;
else
    quant_norm = 0;
end

if(quant_norm)
    [~, sort_perm] = sort(HumanExome.mis_z(J));
    quant_norm_mis_z = sort(HumanChimp.s(I));
    quant_norm_mis_z = quant_norm_mis_z(inv_perm(sort_perm));
    
    figure; plot(HumanExome.mis_z(J), quant_norm_mis_z, '.')
    s_vec_y = normcdf(-HumanExome.mis_z(J)) - 0.1 * HumanExome.mis_z(J);
    s_vec_y = s_vec_y * 2*10^(-4) - 2*10^(-4);
    
    s_vec_x = HumanChimp.s(I);
else % just take selection coefficients
    s_vec_x = HumanChimp.s(I); 
    s_vec_y = -HumanExome.s_het(J);
    s_vec_y_95 = HumanExome.s_het_upper95(J)-HumanExome.s_het_lower95(J);
end

figure; plot(s_vec_x, s_vec_y, '.'); hold on;
[rho, rho_pval] = corr(s_vec_x, s_vec_y);
xlabel('$Human\!-\!Chimp \quad \hat{s}$', 'interpreter', 'latex', 'fontsize', 14);
ylabel('$SFS \quad \hat{s}$', 'interpreter', 'latex', 'fontsize', 14);
[beta, rho2] = polyfit(s_vec_x, s_vec_y, 1);
pos_inds = find(s_vec_x > (-4)*10^(-4));
[beta_pos, rho2_pos] = polyfit(s_vec_x(pos_inds), s_vec_y(pos_inds), 1);
x_vec = min(HumanChimp.s(I)):0.0000001:max(HumanChimp.s(I));
plot(x_vec, x_vec*beta(1)+beta(2), 'r', 'linewidth', 2);
%plot(x_vec, x_vec*beta_pos(1)+beta_pos(2), 'r', 'linewidth', 2);
legend({'data', ['fit (\rho-' num2str(rho, 3) ')']}, 'location', 'northwest', ...
    'fontsize', 14); legend('boxoff');
my_saveas(gcf, fullfile(sfs_figs_dir, 'inter_intra_species_s_hat'), {'epsc', 'pdf', 'jpg'});


% Finally, plot the scatter of the two:
figure;
loglog(-s_vec_x, -s_vec_y, '.'); hold on; % make a scatter plot of s inferred from exome data and from human-chimp comparison
%loglog(min(-s_vec_x):10^(-6):max(-s_vec_x), min(-s_vec_x):10^(-6):max(-s_vec_x), 'r'); % ploy y=x line 
xlabel('$\hat{s}(Ka/Ks)$', 'interpreter', 'latex', 'fontsize', 14); 
ylabel('$\hat{s}_D$', 'interpreter', 'latex', 'fontsize', 14);
% title(['\rho=' num2str(corr(s_vec_x, s_vec_y), 3)]);

special_inds1 = find((-s_vec_x > 0) & ((-s_vec_x < 10^(-5.5)) & (-s_vec_y > 10^(-1)))  ); 
special_inds2 = find((-s_vec_x > 10^(-2)) & (-s_vec_y < 10^(-3)));

loglog(-s_vec_x(special_inds1), -s_vec_y(special_inds1), 'r*');
loglog(-s_vec_x(special_inds2), -s_vec_y(special_inds2), 'g*');
my_saveas(gcf, fullfile(sfs_figs_dir, 'inter_intra_species_s_hat_log'), {'epsc', 'pdf', 'jpg'});






% HERE plot ONLY exomes. Take all (not just intersection!)
s_vec_y = -HumanExome.s_het; % take ALL genes !!
s_vec_y_95 = HumanExome.s_het_upper95-HumanExome.s_het_lower95;


% Plot confidence intervals 
figure; loglog(-s_vec_y, s_vec_y_95 ./ 4, '.'); xlabel('$\hat{s}_{D}$', 'interpreter', 'latex', 'fontsize', 14, 'fontweight', 'bold'); 
ylabel('$\hat{\sigma}(\hat{s}_{D})$', 'interpreter', 'latex', 'fontsize', 14, 'fontweight', 'bold'); 
xlim([10^(-5) 1]); 
my_saveas(gcf, fullfile(sfs_figs_dir, 'PTV_s_hat_vs_std'), {'epsc', 'pdf', 'jpg'});

% plot histogram 
figure; hist(-s_vec_y, 100); xlabel('SFS $\hat{s}_D$', 'interpreter', 'latex'); ylabel('counts'); 

figure; semilogx(sort(-s_vec_y), (1:length(s_vec_y)) ./ length(s_vec_y));  xlabel('SFS $\hat{s}_{D}$', 'interpreter', 'latex'); ylabel('cumulative counts '); % plot cumulative



[h_s, bins_s] = hist(-s_vec_y, 250); 
figure; semilogx(bins_s, smooth(h_s, 5)); xlabel('SFS $\hat{s}_{D}$', 'interpreter', 'latex'); ylabel('density counts '); % plot cumulative
bins_vec = logspace(-5, 0, 100); 
[h_s, bins_s] = hist(-s_vec_y, bins_vec); 
figure; semilogx(bins_s, smooth(h_s, 5)); xlabel('SFS $\hat{s}_{D}$', 'interpreter', 'latex'); ylabel('density counts '); % plot cumulative




