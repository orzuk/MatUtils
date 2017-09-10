% Compute all pairwise correlations between individuals for genotype and phenotype
% Should replicate findings from Yang et al.  Height Nature genetics 2010 paper
%
% Input:
% snp_mat - matrix of snps of all individuals
% trait_vec - vector of trait values for each individual
% trait_type - type of trait (binary or QTL)
%
% Output:
% h_explained - fraction of variance in trait explained liearly by change in genotype
%
function h_explained = visscher_heritability_analysis(snp_mat, trait_vec, trait_type)

genetic_dist_mat = vecs_to_distance(snp_mat); % Compute matrix of pairwise genetic distances
switch trait_type
    case 'QTL'
        trait_dist_mat = vecs_to_distance(trait_vec);
    case 'binary'
        
end



G = genetic_relationship_matrix_internal(snp_mat, snp_type, f_vec); %Build Genetic relashionship matrix

C = corr(genetic_dist_mat(:), trait_dist_mat(:)); % compute correlation coefficient between distances
h_explained = C^2; % fraction of variance explaiend is square of correlation coefficient



