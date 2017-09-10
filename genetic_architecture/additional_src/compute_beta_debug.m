% function compute_beta_debug()

for phenotype_file_name = {'chr1_all_M1_0.05_M2_0.4_h2_0.8.pheno', ...
        'chr1_all_beta_0.1_f_0.4_h2_0.8.txt.pheno', 'chr1_all_lower_0.1_upper_0.2_h2_0.8.pheno'};
    R = ReadDataFile(fullfile('../../common_disease_model/data/fdr/sample_eliana/debug_beta/', ...
        phenotype_file_name{1})); % chr1_all_M1_0.05_M2_0.4_h2_0.8.pheno');
    
    pheno_vec = cell2mat(R.pheno(2:end));
    figure; hold on; hist_density(pheno_vec, 500);
    title(['Histogram of phenotypes. Mean=' num2str(mean(pheno_vec), 3) ', Var=' num2str(var(pheno_vec), 3)]);
    xlabel('phenotypic value'); ylabel('Freq.');
end
