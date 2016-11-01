% Plot comparison of selection coefficients obtained from two different sources
function plot_two_selection_comparison(HumanChimp, HumanExome, sfs_figs_dir)


% Intersect list with the list of selection coefficients from human-exome studies
[common_genes, I, J] = intersect(upper(HumanChimp.Symbol), upper(HumanExome.gene));

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
[rho, rho_pval] = corr(s_vec_x, s_vec_y)
xlabel('$Human-Chimp \quad \hat{s}$', 'interpreter', 'latex', 'fontsize', 14);
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
loglog(min(-s_vec_x):10^(-6):max(-s_vec_x), min(-s_vec_x):10^(-6):max(-s_vec_x), 'r'); % ploy y=x line 
xlabel('Human Chimp s estimator'); ylabel('SFS s estimator');
title(['\rho=' num2str(corr(s_vec_x, s_vec_y), 3)]);

% Plot confidence intervals 
figure; loglog(-s_vec_y, s_vec_y_95, '.'); xlabel('SFS s-hat'); ylabel('SFS-st.d.(s-hat)'); 

% plot histogram 
figure; hist(-s_vec_y, 100); xlabel('SFS s-hat'); ylabel('counts'); 

figure; semilogx(sort(-s_vec_y), 1:length(s_vec_y));  xlabel('SFS s-hat'); ylabel('cumulative counts '); 


