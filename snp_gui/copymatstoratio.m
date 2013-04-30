function [copy_num_vec, allele_ratio_vec] = CopyMatsToRatio(CopyMatA, CopyMatB)

copy_num_vec = CopyMatA+CopyMatB;
allele_ratio_vec = CopyMatA./CopyMatB;

