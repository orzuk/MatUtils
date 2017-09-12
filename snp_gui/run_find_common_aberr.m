%run_find_common_aberr


AssignAllGlobalConstants(); 


chip_type = 'Hind'; 
chip_annot_f = [chip_type '_annot_data_' genome_assembly '.mat'];
genes_db_f = ['refgenes_' genome_assembly '.mat'];
user_dir = '..\data\Leukemia\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load SNP annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNP_annot_struct = load(['../Database/' chip_annot_f]); 
snp_ids = SNP_annot_struct.snp_ids;
chr_loc_vec = SNP_annot_struct.chr_loc_vec;
chr_vec = SNP_annot_struct.chr_vec;
snp_gene_symbols = SNP_annot_struct.snps_gene_symbols;
snps_gene_dist = SNP_annot_struct.snps_gene_dist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load genes description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genes_db_struct = load(['../Database/' genes_db_f]); 
gene_descr = genes_db_struct.gene_descr;
gene_symbols = genes_db_struct.gene_symbols;
% make gene_descr for snp_gene_symbols
snp_gene_descr = snp_gene_symbols_into_descr(snp_gene_symbols, gene_symbols, gene_descr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run function that find stretches and saves output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_names = {'TEL74_5_n';'TEL74_5_d';'TEL76_7_n';'TEL76_7_d'};

%     ;'HD78_9_n';'HD78_9_d';...
%     'TEL80_1_n';'TEL80_1_d';'HD82_3_n';'HD82_3_d';'TEL84_5_n';'TEL84_5_d';'HD86_7_n';'HD86_7_d';...
%     'TEL88_9_n';'TEL88_9_d';'HD90_1_n';'HD90_1_d';'HD92_3_n';'HD92_3_d'}; %;'TEL94_5_n';'TEL94_5_d'};

% deletion: del_amp_flag=1
% amplification: del_amp_flag=2
del_amp_flag = 1; 
thresh_del = 1.6;
thresh_amp = 2.5;
Q = 0.05;

if ~exist(fullfile(user_dir,'output'),'dir')
    mkdir(fullfile(user_dir,'output'));
end
    
find_common_aberr(sample_names, user_dir, chip_type, snp_ids, chr_loc_vec, ...
    chr_vec, snp_gene_symbols, snp_gene_descr, snps_gene_dist, del_amp_flag, ...
    thresh_del, thresh_amp, Q);