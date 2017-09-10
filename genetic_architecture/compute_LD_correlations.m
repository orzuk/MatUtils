% Extract LD correlation matrix from LD file 
%
% Input: 
% LD_file - file with LD information
% 
% Ouptut: 
% LD_corr_mat - matrix of LD values
% 
function LD_corr_mat = compute_LD_correlations( LD_file )

SNPs_data_dir = '../../common_disease_model/data/fdr';
% SNPs_figs_dir = fullfile(SNPs_data_dir, 'figs');

if(~exist('LD_file', 'var') || isempty(LD_file)) % load data
    LD_file = fullfile(SNPs_data_dir, 'plink.ld.chr1.txt'); % file with header
    SNPs_file = fullfile(SNPs_data_dir, 'loci_chr1_with_LD.txt');
end
if(~exist(file_name_to_mat(SNPs_file), 'file'))
    R_SNPs = ReadDataFile(SNPs_file);
    save(file_name_to_mat(SNPs_file), 'R_SNPs');
else
    load(file_name_to_mat(SNPs_file));
end

if(~exist(file_name_to_mat(LD_file), 'file'))
    R_LD = ReadDataFile(LD_file, []);
    save(file_name_to_mat(LD_file), 'R_LD');
else
    R_LD = load(file_name_to_mat(LD_file));
    if(isfield(R_LD, 'R_LD'))
        R_LD = R_LD.R_LD;
    end
end

[~, I, J] = intersect(R_LD.SNP_A, R_SNPs.rsid);
R_LD.f_vec(I) = R_SNPs.all_maf(J); 


%LD_I = cell2vec(R_SNPs.LD_inds);
LD_corr_mat = sparse(III, JJJ, R_LD.R2 .* R_LD.f_vec(III) .* (1-R_LD.f_vec(III)) ./ ...
    (R_LD.f_vec(III) .* (1-R_LD.f_vec(III)))); % store as a sparse matrix 

%num_freq_bins = 10; % how many bins to perform statistics for 
% LD_freq_bins_corr_mat = zeros(num_freq_bins); % global matrix averaging bin frequencies 
