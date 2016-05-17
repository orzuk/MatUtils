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


HumanChimp.w = HumanChimp.Ka ./ HumanChimp.Ks % Compute selection 
bad_inds = union(find(HumanChimp.w > 5), find(isnan(HumanChimp.w)))
HumanChimp.w(bad_inds) = HumanChimp.Ka(bad_inds) ./ HumanChimp.Ki(bad_inds) % Compute selection 
figure; hist(HumanChimp.w, 100); xlabel('w=Ka/Ks'); ylabel('Freq.'); 
num_genes = length(HumanChimp.Ka);

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
plot(S_vec, S_vec + 1, 'r') 

xlabel('S'); ylabel('w = Ka/Ks'); 
%title('Formula from Nielsen and Yang, MBE 2003 for Ka/Ks as function of S'); 
title('Formula: w = S / (1-e^{-S})');  
my_saveas(gcf, fullfile(sfs_figs_dir, 'S_vs_KaKs'), {'epsc', 'pdf', 'jpg'}); 




% Intersect list with the list of selection coefficients from human-exome studies 
[common_genes, I, J] = intersect(HumanChimp.Symbol, SFS.Symbol); 

% Finally, plot the scatter of the two: 
plot(HumanChimp.s(I), SFS.s(J), '.'); % make a scatter plot of s inferred from exome data and from human-chimp comparison
xlabel('Human Chimp s estimator'); ylabel('SFS s estimator'); 



