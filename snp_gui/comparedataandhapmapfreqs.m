% Compare SNPs frequencies from hapmap and from a given data
function [A_data_freq_vec A_hapmap_freq_vec data_snp_ids] = ...
    CompareDataAndHapmapFreqs(user_dir, chip_type, hapmap_dir, hapmap_population)

AssignAllGlobalConstants();

% Get the data snp frequencies
[A_data_freq_vec data_snp_ids] = calc_allele_freq(user_dir, chip_type);


% Get hapmap snp frequencies
load(fullfile(hapmap_dir, pop_str_vec{hapmap_population}, HAPMAP_VERSION, ...
    [chip_type, '_genotypes_chr_' pop_str_vec{hapmap_population} '_' HAPMAP_VERSION '_nr_fwd.mat']));

A_hapmap_freq_vec = -1+zeros(1, length(A_data_freq_vec));

% Run over chromosomes and compare frequencies
for chrom = 1:24
    [snp_intersect I J] = intersect(data_snp_ids, SnpsIDs{chrom});
    A_hapmap_freq_vec(I) = SnpsFreqs{chrom}(J);
end

intersect_inds = find(A_hapmap_freq_vec > -1);

A_data_freq_vec = A_data_freq_vec(intersect_inds);
A_hapmap_freq_vec = A_hapmap_freq_vec(intersect_inds);
data_snp_ids = data_snp_ids(intersect_inds);

figure; plot(A_data_freq_vec, A_hapmap_freq_vec, '.'); xlabel('data freqs'); ylabel('hapmap freqs');



