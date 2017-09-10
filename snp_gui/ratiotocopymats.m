function [CopyMatA CopyMatB] = RatioToCopyMats(copy_num_vec, allele_ratio_vec)

CopyMatB = copy_num_vec ./ (allele_ratio_vec + 1); % We assume sample_ratio is A/B
CopyMatA = CopyMatB .* allele_ratio_vec;
