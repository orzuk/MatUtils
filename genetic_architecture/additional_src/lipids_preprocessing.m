% Analyze lipids data from Eliana.
%
% function lipids_preprocessing()

lipids_file = '../../common_disease_model/data/lipids/lipids_new.txt';
lipids_file2 = '../../common_disease_model/data/lipids/lipids_stage2.txt';
L = ReadDataFile(lipids_file);
L2 = ReadDataFile(lipids_file2);

total_var_vec = [34 15 91].^2; % varianece used to normalize effect sizes (third is just a guess)
traits_vec = {'LDL', 'HDL', 'TG'};
for i=1:length(traits_vec)
    trait = traits_vec{i};  % loop on 3 traits
    trait_inds2 = strfind_cell(L.Secondary, trait);
    trait_inds = strfind_cell(L.Trait, trait)

    [inter_snps I J] = intersect(L.SNP(trait_inds2), L2.MarkerName);
    trait_beta2 = zeros(1, length(I)); 
    trait_beta2(I) = L2.EffectSize(J); 
    trait_beta = [L.Effectsize(trait_inds)' trait_beta2]';
    trait_MAF = L.MAF([trait_inds trait_inds2]);
    trait_inds = [trait_inds trait_inds2]
    
    trait_V = beta_to_variance_explained(trait_beta, trait_MAF, total_var_vec(i), 'diploid');
    trait_beta = trait_beta ./ sqrt(total_var_vec(i));    
    trait_SNP = L.SNP(trait_inds);
    trait_chr = L.Chr(trait_inds);
    trait_locus = L.Locus(trait_inds);
    trait_pval = L.p_value(trait_inds); % wrong for secondary traits but what the hell .. 
    R = [trait_chr trait_SNP trait_locus trait_pval num2cell(trait_MAF)  num2cell(trait_beta)];
    R = [{'chr', 'snp-id', 'gene', 'p-val', 'MAF', 'beta'}' R']'
    savecellfile(R, [lipids_file(1:end-4) '_' trait '.txt']);
end



