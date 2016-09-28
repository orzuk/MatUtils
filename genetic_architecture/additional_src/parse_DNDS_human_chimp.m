% Parse dn-ds data from .. then compare to selection coefficient s
% inferred from exon sequencing data 

AssignGeneralConstants
AssignRVASConstants;
DNDS_file_name = 'C:\research\RVAS\Data\nature04072-s6.txt'

% Read human-chimp file:
if(exist(file_name_to_mat(DNDS_file_name), 'file'))
    HumanChimp = load(file_name_to_mat(DNDS_file_name));
else
    [HumanChimp, R] = ReadDataFile(DNDS_file_name, [], 1, [], tab);
end
HumanChimp.Symbol = num2str_cell(HumanChimp.Symbol); % correct problems in gene names (position appears instead)
num_genes = length(HumanChimp.Ka);

HumanChimp.w = HumanChimp.Ka ./ HumanChimp.Ks % Compute selection 
bad_inds = union(find(HumanChimp.w > 5), find(isnan(HumanChimp.w)))
good_inds = setdiff(1:num_genes, bad_inds);

HumanChimp.w(bad_inds) = HumanChimp.Ka(bad_inds) ./ HumanChimp.Ki(bad_inds) % Compute selection 
figure; hist(HumanChimp.w, 100); xlabel('w=Ka/Ks'); ylabel('Freq.'); 

% fit s by solving : Ka/Ks = S / (1-e^{-S})   from Nielsen&Yang
HumanChimp.s = zeros(num_genes,1); 

for i=1:num_genes
    if(mod(i, 100) == 0)
        fit_s_i = i
    end
    HumanChimp.s(i) = fzero(@(S) S - (1-exp(-S))*max(10^(-6),HumanChimp.w(i)), 2*log(max(10^(-6),HumanChimp.w(i)))); % solve with close-by initial condition 
end
HumanChimp.s = HumanChimp.s ./ (4*N_eff); 

figure; hist(HumanChimp.s, 100); xlabel('s-hat (Human-Chimp)'); ylabel('Freq.'); 
my_saveas(gcf, fullfile(sfs_figs_dir, 'Inter_vs_Intra_species_selection', 'S_hat_distribution_human_chimp'), {'epsc', 'pdf', 'jpg'}); 



% Plot theoretical Curve from Nielsen&Yang
S_vec = -5:0.01:5; w_vec = S_vec ./ (1-exp(-S_vec));
figure; hold on; plot(S_vec, w_vec); 
xlabel('S'); ylabel('w = Ka/Ks'); 
%title('Formula from Nielsen and Yang, MBE 2003 for Ka/Ks as function of S'); 
title('Formula: w = S / (1-e^{-S})');  
my_saveas(gcf, fullfile(sfs_figs_dir, 'S_vs_KaKs'), {'epsc', 'pdf', 'jpg'}); 



% Get SFS selection scores 
use_ExAC_constraint=1; 
if(use_ExAC_constraint) % use constraint from Samocha et al. 
     Constraint_file_name = 'C:\research\RVAS\Data\SiteFrequencySpectra\ExAC\ConstraintSamocha\fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'; 
    % Read constraint file:
    if(exist(file_name_to_mat(Constraint_file_name), 'file'))
        HumanExome = load(file_name_to_mat(Constraint_file_name));
    else
        [HumanExome, R] = ReadDataFile(Constraint_file_name, [], 1, [], tab);
    end
    HumanExome.num_genes = length(HumanExome.gene); 
    figure; hold on; 
    plot(sort(HumanExome.syn_z), (1:HumanExome.num_genes)./HumanExome.num_genes); 
    plot(sort(HumanExome.mis_z), (1:HumanExome.num_genes)./HumanExome.num_genes, 'r'); 
    plot(sort(HumanExome.lof_z), (1:HumanExome.num_genes)./HumanExome.num_genes, 'g'); 
    legend({'Synonymous', 'Missense', 'LOF'}); legend('boxoff'); xlabel('Z-score'); ylabel('Cumulative freq.'); 
     
     
else % use our s-hat estimator 
    
    
end

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



