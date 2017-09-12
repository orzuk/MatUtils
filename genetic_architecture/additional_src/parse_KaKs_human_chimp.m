% Parse DN-DS data from Chimp-genome paper (Nature).
% Then compare to selection coefficient s inferred from exon sequencing data

AssignGeneralConstants
AssignRVASConstants;

DNDS_file_name = 'C:\research\RVAS\Data\nature04072-s6.txt'; % file with human-chimp DN, DS values
% Get SFS selection scores
exome_constraint='ExAC_PTV'; % 'ExAC_PTV'
switch exome_constraint
    case 'ExAC_Samocha'  % use constraint from Samocha et al.
        Constraint_file_name = 'C:\research\RVAS\Data\SiteFrequencySpectra\ExAC\ConstraintSamocha\fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt';
    case 'ExAC_PTV' % protein truncation varinats from Shamil's biorxiv
        Constraint_file_name = 'C:\research\RVAS\Data\SiteFrequencySpectra\ExAC\PTV\ProteinTruncatingVarinats_selection_075523-1.txt'; % filled file name!!!
    case 'ExAC_TwoClass'; % s_hat estimator from our model using missense+stops
        Constraint_file_name = 'C:\research\RVAS\Data\SiteFrequencySpectra\ExAC\out\GeneByGene.txt'; % fill file name!!! (for each population!)
end


% Read human-chimp file:
if(exist(file_name_to_mat(DNDS_file_name), 'file'))
    HumanChimp = load(file_name_to_mat(DNDS_file_name));
else
    [HumanChimp, R] = ReadDataFile(DNDS_file_name, [], 1, [], tab);
end
HumanChimp.Symbol = num2str_cell(HumanChimp.Symbol); % correct problems in gene names (position appears instead)

%if(~isfield(HumanChimp, 's')) % compute selection coefficient (once)
    [HumanChimp.s, HumanChimp.w] = KaKs_to_selection_coefficient(HumanChimp.Ka, HumanChimp.Ks); % Compute selection
    num_genes = length(HumanChimp.Ka);
    HumanChimp.bad_inds = union(find(HumanChimp.w > 5), find(isnan(HumanChimp.w)))
    HumanChimp.good_inds = setdiff(1:num_genes, HumanChimp.bad_inds);
    if(~isempty(HumanChimp.bad_inds))
        [HumanChimp.s(HumanChimp.bad_inds), HumanChimp.w(HumanChimp.bad_inds)] = ...
            KaKs_to_selection_coefficient(HumanChimp.Ka(HumanChimp.bad_inds), HumanChimp.Ki(HumanChimp.bad_inds));  % Compute selection
    end
    save(file_name_to_mat(DNDS_file_name), '-struct', 'HumanChimp'); % save again, with selected coefficients
%end

figure; hist(HumanChimp.w, 100); xlabel('w=Ka/Ks'); ylabel('Freq.');
figure; hist(HumanChimp.s, 100); xlabel('s-hat (Human-Chimp)'); ylabel('Freq.');
my_saveas(gcf, fullfile(sfs_figs_dir, 'Inter_vs_Intra_species_selection', 'S_hat_distribution_human_chimp'), {'epsc', 'pdf', 'jpg'});
% Plot theoretical Curve from [Nielsen&Yang]
S_vec = -5:0.01:5; w_vec = S_vec ./ (1-exp(-S_vec));
figure; hold on; plot(S_vec, w_vec);
xlabel('S'); ylabel('w = Ka/Ks');
title('Formula: w = S / (1-e^{-S})');   %title('Formula from Nielsen and Yang, MBE 2003 for Ka/Ks as function of S');
my_saveas(gcf, fullfile(sfs_figs_dir, 'S_vs_KaKs'), {'epsc', 'pdf', 'jpg'});


% Read constraint file:
if(exist(file_name_to_mat(Constraint_file_name), 'file'))
    HumanExome = load(file_name_to_mat(Constraint_file_name));
else
    [HumanExome, R] = ReadDataFile(Constraint_file_name, [], 1, [], tab);
end
if(isfield(HumanExome, 'gene'))
    HumanExome.num_genes = length(HumanExome.gene);
else
    HumanExome.num_genes = length(HumanExome.gene_symbol);
    HumanExome.gene = HumanExome.gene_symbol;
end

if(isfield(HumanExome, 'syn_z')) % plot SFS
    figure; hold on;
    plot(sort(HumanExome.syn_z), (1:HumanExome.num_genes)./HumanExome.num_genes);
    plot(sort(HumanExome.mis_z), (1:HumanExome.num_genes)./HumanExome.num_genes, 'r');
    plot(sort(HumanExome.lof_z), (1:HumanExome.num_genes)./HumanExome.num_genes, 'g');
    legend({'Synonymous', 'Missense', 'LOF'}); legend('boxoff'); xlabel('Z-score'); ylabel('Cumulative freq.');
end

% Perform comparison of selection coefficients inferred from different sources
plot_two_selection_comparison(HumanChimp, HumanExome, sfs_figs_dir);


