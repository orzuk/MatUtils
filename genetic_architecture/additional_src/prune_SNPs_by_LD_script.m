% Prune table of SNP pvals based on LD structure
%
% Input:
% SNPs_file - file with data on all SNPs
% LD_file - file with pairwise LD data
% GWAS_results_file - file with p-values, FDR etc. for each SNP
% causal_SNPs_file - New: add a possibility of the true causal SNPs. This is used to remove stuff in LD with them
% operation_flag - 1: parse LD files, 0 - prune
% save_txt_flag - 1: save data also in .txt tab-delimited format (default, slow). 0 - save only in .mat format
% trait_type - quantitative (default) or binary (disease)
% perform_deconvolution - flag saying if to do deconvolution or not (default, NOT)
%
function prune_SNPs_by_LD_script(SNPs_file, LD_file, GWAS_results_file, causal_SNPs_file, ...
    operation_flag, save_txt_flag, trait_type, perform_deconvolution)

AssignGeneralConstants;
% Never add path
% if isdeployed
%     SetPathScript;
% end

if(nargin == 0) % no inputs. Display parameters
    sprintf(['Usage: ./prune_SNPs_by_LD_script_exe SNPs_file LD_file GWAS_results_file causal_SNPs_file operation_flag\n\n' ...
        'Parameters:\n' ...
        'SNPs_file - large file with information on SNPs (needed only for pre-processing or to get allele frequencies)\n' ...
        'LD_file - file containing all pairwise-LD information in a convenient form\n' ...
        'GWAS_results_file - file containing GWAS results, with effect size, p-value and fdr\n' ...
        'causal_SNPs_file -  file containing true causal SNPs and their effect size (when simulating data)\n' ...
        'operation_flag - what to do: 0. prune_snps (default), 1. preprocess_ld 2. prepare LD matrices \n\n' ...
        'Output pruned file (for the default prune_snps operation) will be in the file: [GWAS_results-without-suffix ''_with_LD_pruned.txt'']\n\n' ...
        'Example: ./prune_SNPs_by_LD_script_exe ''''  ' ...
        '/seq/orzuk/common_disease_model/data/fdr/PLINK/CEU_chr1.ld_preprocessed.txt ' ...
        '/seq/orzuk/common_disease_model/data/fdr/chr1_20_causal_snps_bin_40_50.txt ' ...
        '/seq/orzuk/common_disease_model/data/fdr/causal_loci_chr1.txt prune_snps\n\n' ...
        'For this example, output will be in: /seq/orzuk/common_disease_model/data/fdr/chr1_20_causal_snps_bin_40_50_with_LD_pruned.txt'])
    return;
end

%%% We don't need to load the causal files here.
if(~isempty(causal_SNPs_file)) % see if we have a causal file name
    if(~exist(file_name_to_mat(causal_SNPs_file), 'file'))
        C = loadcellfile(causal_SNPs_file);
        if(isempty(C))
            return;
        end
        %    save(file_name_to_mat(causal_SNPs_file), 'C'); % save in .mat format
    else % just load .mat file
        load(file_name_to_mat(causal_SNPs_file));
    end
    if(exist('C', 'var') & (isempty(C))) % indicative of empty
        return;
    end
end % if we have a causal-SNPs file

switch machine
    case PC
        SNPs_data_dir = '../../common_disease_model/data/fdr';
    case UNIX
        SNPs_data_dir = '/seq/orzuk/common_disease_model/data/fdr';
end
run_str = 'all'; % 'all'; % 'small'; % 'all';
% LD_thresh = 0.5;

if(~exist('operation_flag', 'var') || isempty(operation_flag)) % set default: compare LD
    operation_flag = 0;
end
if(~exist('save_txt_flag', 'var') || isempty(save_txt_flag)) % set default: compare LD
    save_txt_flag = 1;
end
if(~exist('perform_deconvolution', 'var') || isempty(perform_deconvolution))
    perform_deconvolution = 1;
end

% if(exist('SNPs_file', 'var'))
%     input_SNPs_file = SNPs_file
% end
if(~exist('SNPs_file', 'var') || isempty(SNPs_file))
    switch run_str
        case 'all'  % run the full data (chromosome 1)
            SNPs_file = fullfile(SNPs_data_dir, 'loci_chr1.txt');
        case 'small' % small test run
            SNPs_file = fullfile(SNPs_data_dir, 'loci_small.txt');
    end
end
if(~exist('LD_file', 'var') || isempty(LD_file))
    switch run_str
        case 'all'  % run the full data (chromosome 1)
            LD_file = fullfile(SNPs_data_dir, 'plink.ld.chr1.txt'); % file with header
        case 'small' % small test run
            LD_file = fullfile(SNPs_data_dir, 'plink.ld.small.txt');
    end
    
end
if(~exist('GWAS_results_file', 'var') || isempty(GWAS_results_file))
    GWAS_results_file = fullfile(SNPs_data_dir, 'chr1_1.snptest.out.wfdr.foror.txt');
end

switch operation_flag % What operation should be doing
    case {0, 'prune_snps'} % prune LD
        prune_snps(SNPs_file, LD_file, GWAS_results_file, causal_SNPs_file, SNPs_data_dir, ...
            trait_type, save_txt_flag, perform_deconvolution);
        
    case {1, 'preprocess_ld'} % parse LD files. Do not need the SNP files
        sprintf('preprocess LD')
        preprocess_LD(LD_file, save_txt_flag);
        
    case {2, 'prepare_LD_matrices'}  % prepare deconvolution matrices
        prepare_LD_matrix(LD_file);
end % switch operation_flag



% Internal function: Plot var explained



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal function: prune snps according to LD
%
% Input:
% SNPs_file - (optional) file containing allele-frequency information for SNPs
% LD_file - file with pairwise ld information (alreadey pre-processed)
% GWAS_results_file - file with effect size and p-value for each snp
% causal_SNPs_file - file with real effect sizes for causal SNPs
% trait_type - quantitative (default) or disease
% save_txt_flag - flag saying if yo save output also in .txt format
% perform_deconvolution - flag saying if to perform deconvolution or NOT (default)
%
function prune_snps(SNPs_file, LD_file, GWAS_results_file, causal_SNPs_file, SNPs_data_dir, ...
    trait_type, save_txt_flag, perform_deconvolution)

% perform_deconvolution=0; % dont try deconvolution for now
SNPs_figs_dir = fullfile(SNPs_data_dir, 'figs');

run_str = dir_from_file_name(GWAS_results_file); run_str = remove_dir_from_file_name(run_str(1:end-1));

if(exist(file_name_to_mat(LD_file), 'file'))  % load R_LD
    R_LD = load(file_name_to_mat(LD_file), 'rsid', 'LD_inds'); % , 'R2'); % no need for R2
    if(isempty(fields(R_LD)))
        R_LD = load(file_name_to_mat(LD_file));
    end
    if(isfield(R_LD, 'R_LD'))
        R_LD = R_LD.R_LD;
    end
else
    R_LD = ReadDataFile(LD_file);
end
%        load([remove_suffix_from_file_name(SNPs_file) '_with_LD.mat']); % load R_SNPs
if(exist(file_name_to_mat(GWAS_results_file), 'file'))
    GWAS_RES = load(file_name_to_mat(GWAS_results_file), 'rsid', 'beta', 'p', 'pval', 'fdr', 'FDR', 'f_vec'); % load only needed fields to save time
    if(isempty(fields(GWAS_RES)))
        GWAS_RES = load(file_name_to_mat(GWAS_results_file)); % load everything
        if(isfield(GWAS_RES, 'GWAS_RES'))
            GWAS_RES = GWAS_RES.GWAS_RES;
        end
    end
else
    reading_file = GWAS_results_file
    reading_output_file = file_name_to_mat(GWAS_results_file);
    fid = fopen(GWAS_results_file, 'r');  % deal with possibly empty first line
    tmp_line = fgetl(fid);
    if(isempty(tmp_line))
        skip_lines = 1;
    else
        skip_lines = 0;
    end
    GWAS_RES = ReadDataFile(GWAS_results_file, file_name_to_mat(GWAS_results_file), [], skip_lines);
    if(isempty(GWAS_RES)) % empty file - nothing to do
        return;
    end
    save(file_name_to_mat(GWAS_results_file), '-struct', 'GWAS_RES');
end
if(~isfield(GWAS_RES, 'FDR')) % set fdr field name
    GWAS_RES.FDR = GWAS_RES.fdr;
end
GWAS_RES = my_rmfield(GWAS_RES, 'fdr');
if(~isfield(GWAS_RES, 'pval')) % set p-value field name
    GWAS_RES.pval = GWAS_RES.p;
end
GWAS_RES = my_rmfield(GWAS_RES, 'p');
GWAS_RES = sort_struct_by_field(GWAS_RES, 'FDR', 'descend');
num_snps = length(GWAS_RES.rsid);
[intersect_SNPs I J] = intersect(R_LD.rsid, GWAS_RES.rsid);

if(~isfield(GWAS_RES, 'f_vec') || (max(GWAS_RES.f_vec) < 0)) % fill allele frequencies
    GWAS_RES.f_vec = zeros(length(GWAS_RES.rsid), 1);
    if(isfield(R_LD, 'f_vec'))
        GWAS_RES.f_vec(J) = R_LD.f_vec(I);
    else % need to load also SNPs file
        %            now_SNPs_file_is = SNPs_file
        if(~exist(file_name_to_mat(SNPs_file), 'file'))  % load R_SNPs
            fid = fopen(SNPs_file, 'r');  % deal with possibly empty first line
            tmp_line = fgetl(fid);
            if(isempty(tmp_line))
                skip_lines = 1;
            else
                skip_lines = 0;
            end
            R_SNPs = ReadDataFile(SNPs_file, [], [], skip_lines);
            save(file_name_to_mat(SNPs_file), 'R_SNPs');
        else
            load(file_name_to_mat(SNPs_file));
        end
        [~, I_SNPs J_SNPs] = intersect(R_SNPs.rsid, GWAS_RES.rsid);
        if(~isfield(R_SNPs, 'f_vec'))
            R_SNPs.f_vec = R_SNPs.all_maf;
        end
        GWAS_RES.f_vec(J_SNPs) = R_SNPs.f_vec(I_SNPs);
    end
end % if isfield f_vec
GWAS_RES.f_vec = min(GWAS_RES.f_vec, 1 - GWAS_RES.f_vec); % transform to minor-allele-frequencies .. 

if(exist('allele_freq_file', 'var')) % Load allele frequencies from hapmap samples
    if(exist(allele_freq_file, 'file'))
        
    end
end

J_inv = zeros(num_snps,1);
for i=1:length(J)
    J_inv(J(i))=i;
end
%        J_inv = inv_perm(J); I_inv = inv_perm(I);
GWAS_RES.keep_vec = ones(num_snps,1);
for i=1:num_snps % length(J) % Loop on SNPs according to p-value
    if(mod(i,1000) == 0)
        prune_snp_i = i
    end
    if(J_inv(i) > 0)
        cur_ind = I(J_inv(i));
        if(GWAS_RES.keep_vec(i) == 1) % didn't exclue this one yet
            GWAS_RES.keep_vec(R_LD.LD_inds{cur_ind}) = 0;
        end
    end
end
J_comp = setdiff(1:num_snps,J); GWAS_RES.keep_vec(J_comp)=0; % get rid of snps without p-vals
GWAS_RES.keep_vec = GWAS_RES.keep_vec(1:num_snps);

keep_inds_GWAS_RES = find(GWAS_RES.keep_vec);
%         [~, keep_inds_GWAS_RES] = intersect(GWAS_RES.rsid, R_SNPs.rsid(keep_SNPs_inds_vec));
%         R_SNPs_pruned = struct_by_inds(R_SNPs, keep_SNPs_inds_vec);
%         GWAS_RES.keep_vec = zeros(length(GWAS_RES.rsid),1); GWAS_RES.keep_vec(keep_inds_GWAS_RES) = 1;
GWAS_RES_pruned = struct_by_inds(GWAS_RES, keep_inds_GWAS_RES);
GWAS_RES_pruned = sort_struct_by_field(GWAS_RES_pruned, 'FDR', 'descend');
save([remove_suffix_from_file_name(GWAS_results_file) '_with_LD_pruned.mat'], ...
    'GWAS_RES', 'GWAS_RES_pruned'); % save in binary .mat format, no we have pruning information
if(save_txt_flag) % save as .txt
    WriteDataFile(GWAS_RES, ...
        [remove_suffix_from_file_name(GWAS_results_file) '_with_LD_pruned.txt']); % save in .txt format
end

if(~isfield(GWAS_RES, 'n_samples'))
    GWAS_RES.n_samples = str2nums(GWAS_results_file);
    if(~isempty(GWAS_RES.n_samples))
        GWAS_RES.n_samples = GWAS_RES.n_samples(end); % convention: last word in file name is # individuals
    else
        GWAS_RES.n_samples = 2000; % set sample size
    end
end
n_samples = GWAS_RES.n_samples;

switch trait_type
    case 'quantitative'
        GWAS_RES.var_explained = beta_to_heritability(GWAS_RES.beta, ...
            GWAS_RES.f_vec, 1, 'diploid'); %  n_samples); % correct for finite sample size bias in var. explained
    case 'disease' % need to transfer to beta
        if(~exist('prevalence', 'var'))
            prevalence = 0.01; % set default
        end
        [~, GWAS_RES.beta] = genetic_relative_risk_to_beta(GWAS_RES.f_vec, GWAS_RES.beta, prevalence);
        GWAS_RES.var_explained = beta_to_heritability(GWAS_RES.beta, ...
            GWAS_RES.f_vec, 1, 'diploid'); %  n_samples); % correct for finite sample size bias in var. explained
end


GWAS_RES.var_explained_times_fdr = GWAS_RES.var_explained .* ...
    GWAS_RES.FDR;
%         GWAS_RES.var_explained2 = beta_to_heritability(GWAS_RES.pheno_frequentist_add_score_beta_1, ...
%             GWAS_RES.f_vec, 1, 'diploid', n_samples)  .* GWAS_RES.pheno_frequentist_add_score_local_fdr_1; %  n_samples);


%         2.*GWAS_RES.pheno_frequentist_add_score_beta_1 .^ 2 .* ...
%             GWAS_RES.f_vec .* (1-GWAS_RES.f_vec) .* GWAS_RES.pheno_frequentist_add_score_local_fdr_1;

small_pvals_inds = intersect(find(GWAS_RES.pval < 0.5), ...
    find(min(GWAS_RES.f_vec, 1-GWAS_RES.f_vec) > 0.01)); % set a very liberal threshold, also remove low MAF
keep_inds_GWAS_RES = intersect(keep_inds_GWAS_RES, small_pvals_inds);
GWAS_RES.f_vec = min(GWAS_RES.f_vec, 1-GWAS_RES.f_vec); % get MAF
res = 1/50;
num_bins = [0:res:0.5]; num_bins(1)=num_bins(1)-eps; num_bins(end)=num_bins(end)+eps;% acutal number is plus 1
[all_var_hist all_var_bins] = weighted_hist(GWAS_RES.f_vec(small_pvals_inds), ...
    GWAS_RES.var_explained_times_fdr(small_pvals_inds), num_bins);
[pruned_var_hist pruned_var_bins] = weighted_hist(GWAS_RES.f_vec(keep_inds_GWAS_RES), ...
    GWAS_RES.var_explained_times_fdr(keep_inds_GWAS_RES), num_bins);
V_all = sum(GWAS_RES.var_explained_times_fdr(small_pvals_inds));
V_pruned = sum(GWAS_RES.var_explained_times_fdr(keep_inds_GWAS_RES));
model_str = remove_suffix_from_file_name(remove_dir_from_file_name(GWAS_results_file));

if(perform_deconvolution) % New: perform deconvolution to remove LD effects
    GWAS_RES.beta_deconvolved = LD_deconvolution(GWAS_RES.beta, R_LD); % perform deconvulution
    GWAS_RES.pval_deconvolved = XXX; % compute p-value to the deconvolved signal
    GWAS_RES.loc_fdr_deconvolved = run_fdr(GWAS_RES.beta_deconvolved);  % Compute new fdr on the deconvolved signal !!!
end % deconvulution

if ~isdeployed
    all_var_bin_size = all_var_bins(2)-all_var_bins(1);
    figure; hold on; bar(all_var_bins, all_var_hist ./ sum(all_var_hist), 0.2);
    bar(pruned_var_bins+all_var_bin_size/4, pruned_var_hist ./ sum(pruned_var_hist), 0.2, 'r');
    xlabel('MAF'); ylabel('Var. Explained');
    legend({['All (V=' num2str(V_all, 4) ')'], ['Pruned (V=' num2str(V_pruned, 4) ')']}, 2);
    my_title(['Estimated Variance Explained, ' model_str]);
    all_y_lim = get(gca, 'ylim')
    my_saveas(gcf, fullfile(SNPs_figs_dir, [remove_dir_from_file_name(GWAS_results_file) '_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
end
if(exist(causal_SNPs_file, 'file') && (~isempty(causal_SNPs_file))) % seperate into stuff in LD and not from causal SNPs
    fid = fopen(causal_SNPs_file, 'r');  % deal with possibly empty first line
    tmp_line = fgetl(fid);
    if(isempty(tmp_line))
        skip_lines = 1;
    else
        skip_lines = 0;
    end
    if(~exist(file_name_to_mat(causal_SNPs_file), 'file'))
        R_causal = ReadDataFile(causal_SNPs_file, [], [], skip_lines); % load causal snps
    else
        R_causal = load(file_name_to_mat(causal_SNPs_file));
    end
    if(~isfield(R_causal, 'f_vec'))  % rename
        R_causal.f_vec = R_causal.maf;
    end
    R_causal.f_vec = min(R_causal.f_vec, 1-R_causal.f_vec); % tranform to minor-allele-frequencies 
    num_causal_snps = length(R_causal.rsid)
    [causal_SNPs I_causal J_causal] = intersect(R_LD.rsid, R_causal.rsid); % get indices of causal SNPs
    [~, I_causal2, J_causal2] = intersect(R_causal.rsid, GWAS_RES.rsid);
    GWAS_RES.f_vec(J_causal2) = R_causal.f_vec(I_causal2); % change allele frequencies of causal SNPs
    causal_LD_SNP_inds = cell2vec(str2nums_cell(R_LD.LD_inds(I_causal)));
    if(~isfield(R_causal, 'beta')) % set effect sizes as equal
        R_causal.beta = repmat( sqrt(0.8 / (num_causal_snps * 0.5) ), num_causal_snps, 1);
    end
    R_causal.var_explained = 2 * R_causal.beta.^ 2 .* R_causal.f_vec .* (1-R_causal.f_vec);
    [only_causal_var_hist only_causal_var_bins] = weighted_hist(GWAS_RES.f_vec(J_causal2), ...
        R_causal.var_explained(I_causal2) .* GWAS_RES.FDR(J_causal2), num_bins);
    V_true = sum(R_causal.var_explained);
    V_true_FDR = sum(R_causal.var_explained(I_causal2) .* GWAS_RES.FDR(J_causal2));
    only_causal_bin_size = only_causal_var_bins(2)-only_causal_var_bins(1);
    model_str = strdiff(model_str, '_wheader');
    if ~isdeployed
        figure; hold on;
        subplot(2,2,1); hold on;
        bar(all_var_bins, all_var_hist ./ sum(all_var_hist), 0.2);
        bar(pruned_var_bins+all_var_bin_size/4, pruned_var_hist ./ sum(pruned_var_hist), 0.2, 'r');
        xlabel('MAF'); ylabel('Var. Explained');
%        legend({['All (V=' num2str(V_all, 4) ')'], ['Pruned (V=' num2str(V_pruned, 4) ')']}, 2);
        my_title(['Estimated Var. Expl., ' model_str]);
        bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist./sum(only_causal_var_hist), 0.2, 'g'); % plot causal snps histogram
        all_y_lim = get(gca, 'ylim')
        legend({['All (V=' num2str(V_all, 3) ')'], ...
            ['Pruned (V=' num2str(V_pruned, 3) ')'], ...
            ['True (V=' num2str(V_true_FDR, 3) ' of ' num2str(V_true, 3) ')']}, 1, 'fontsize', 6);
        xlim([0 0.51]);
        %         my_saveas(gcf, fullfile(SNPs_figs_dir, ...
        %             [remove_dir_from_file_name(GWAS_results_file) '_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
    end
    
    [~, ~, causal_LD_SNP_inds2] = intersect(R_LD.rsid(causal_LD_SNP_inds), GWAS_RES.rsid);
    causal_LD_SNP_inds2 = union(causal_LD_SNP_inds2, J_causal2);
    non_causal_LD_SNP_inds2 = setdiff(1:length(GWAS_RES.rsid), causal_LD_SNP_inds2);
    non_causal_LD_SNP_inds2 = intersect(non_causal_LD_SNP_inds2, small_pvals_inds);
    causal_LD_SNP_inds2 = intersect(causal_LD_SNP_inds2, small_pvals_inds);
    causal_keep_inds_GWAS_RES = intersect(keep_inds_GWAS_RES, causal_LD_SNP_inds2);
    non_causal_keep_inds_GWAS_RES = intersect(keep_inds_GWAS_RES, non_causal_LD_SNP_inds2);
    non_causal_keep_inds_GWAS_RES_vec = zeros(num_snps,1); non_causal_keep_inds_GWAS_RES_vec(non_causal_keep_inds_GWAS_RES)=1;
    
    [sorted_non_causal_pvals sort_non_causal_perm] =  sort(GWAS_RES.pval);
    var_explained_times_fdr_cumsum = cumsum( ...
        GWAS_RES.var_explained_times_fdr(sort_non_causal_perm) .* non_causal_keep_inds_GWAS_RES_vec);
    total_var_explained_non_causal_allowed = 0.10;
    pval_cutoff = find(var_explained_times_fdr_cumsum >= total_var_explained_non_causal_allowed, 1);
    pval_cutoff = sorted_non_causal_pvals(pval_cutoff);
    if(isempty(pval_cutoff)) % everything was below cutoff
        pval_cutoff = 0.5;
    end
    small_pvals_inds = intersect(small_pvals_inds, find(GWAS_RES.pval < pval_cutoff)); % set a strict threshold
    keep_inds_GWAS_RES = find(GWAS_RES.keep_vec);
    keep_inds_GWAS_RES = intersect(keep_inds_GWAS_RES, small_pvals_inds);
    [all_var_hist all_var_bins] = weighted_hist(GWAS_RES.f_vec(small_pvals_inds), ...
        GWAS_RES.var_explained_times_fdr(small_pvals_inds), num_bins);
    [pruned_var_hist pruned_var_bins] = weighted_hist(GWAS_RES.f_vec(keep_inds_GWAS_RES), ...
        GWAS_RES.var_explained_times_fdr(keep_inds_GWAS_RES), num_bins);
    V_all = sum(GWAS_RES.var_explained_times_fdr(small_pvals_inds));
    V_pruned = sum(GWAS_RES.var_explained_times_fdr(keep_inds_GWAS_RES));
    model_str = remove_suffix_from_file_name(remove_dir_from_file_name(GWAS_results_file));
    if ~isdeployed
        subplot(2,2,2);  % Change ! put 3 figures in the same file
        hold on;
        bar(all_var_bins, all_var_hist./sum(all_var_hist), 0.2);
        bar(pruned_var_bins+all_var_bin_size/4, pruned_var_hist./sum(pruned_var_hist), 0.2, 'r');
        xlabel('MAF'); ylabel('Var. Explained');
%        legend({['All (V=' num2str(V_all, 4) ')'], ['Pruned (V=' num2str(V_pruned, 4) ')']}, 2);
        my_title(['Only p-val < ' num2str(pval_cutoff, 3)]);
        bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist./sum(only_causal_var_hist), 0.2, 'g');
        legend({['All (V=' num2str(V_all, 3) ')'], ...
            ['Pruned (V=' num2str(V_pruned, 3) ')'], ...
            ['True (V=' num2str(V_true_FDR, 3) ' of ' num2str(V_true, 3) ')']}, 1, 'fontsize', 6);
        xlim([0 0.51]);
        ylim(all_y_lim); % set y lim as in previous fig. (why would we do that?)
        %         my_saveas(gcf, fullfile(SNPs_figs_dir, ...
        %             [remove_dir_from_file_name(GWAS_results_file) '_var_explained_with_pruning_high_res_strict_pval']), {'epsc', 'jpg', 'pdf'});
    end
    
    
    [causal_var_hist causal_var_bins] = weighted_hist(GWAS_RES.f_vec(causal_LD_SNP_inds2), ...
        GWAS_RES.var_explained_times_fdr(causal_LD_SNP_inds2), num_bins);
    [non_causal_var_hist non_causal_var_bins] = weighted_hist(GWAS_RES.f_vec(non_causal_LD_SNP_inds2), ...
        GWAS_RES.var_explained_times_fdr(non_causal_LD_SNP_inds2), num_bins);
    [pruned_causal_var_hist pruned_causal_var_bins] = weighted_hist(GWAS_RES.f_vec(causal_keep_inds_GWAS_RES), ...
        GWAS_RES.var_explained_times_fdr(causal_keep_inds_GWAS_RES), num_bins);
    [pruned_non_causal_var_hist pruned_non_causal_var_bins] = weighted_hist(GWAS_RES.f_vec(non_causal_keep_inds_GWAS_RES), ...
        GWAS_RES.var_explained_times_fdr(non_causal_keep_inds_GWAS_RES), num_bins);
    V_causal = sum(GWAS_RES.var_explained_times_fdr(causal_LD_SNP_inds2));
    V_non_causal = sum(GWAS_RES.var_explained_times_fdr(non_causal_LD_SNP_inds2));
    V_pruned_causal = sum(GWAS_RES.var_explained_times_fdr(causal_keep_inds_GWAS_RES));
    V_pruned_non_causal = sum(GWAS_RES.var_explained_times_fdr(non_causal_keep_inds_GWAS_RES));
    
    if ~isdeployed
        plot_by_LD = 0;
        if(plot_by_LD)
            figure;   % New: Put all in the same file   % figure;
            hold on; bar(causal_var_bins, causal_var_hist./sum(causal_var_hist), 0.2);
            bar(pruned_causal_var_bins+all_var_bin_size/4, pruned_causal_var_hist./sum(pruned_causal_var_hist), 0.2, 'r');
            bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist./sum(only_causal_var_hist), 0.2, 'g');
            xlabel('MAF'); ylabel('Var. Explained');
            legend({['All LD-causal (V=' num2str(V_causal,4) ')'], ...
                ['Pruned LD-causal (V=' num2str(V_pruned_causal,4) ')'], ...
                ['True (V=' num2str(V_true_FDR, 4) ' of ' num2str(V_true, 4) ')']}, 2, 'fontsize', 8);
            my_title(['Estimated Variance Explained, ' model_str]);
            xlim([0 0.51]); ylim(all_y_lim);
            %         my_saveas(gcf, fullfile(SNPs_figs_dir, [remove_dir_from_file_name(GWAS_results_file) ...
            %             '_allele_freq_var_dist']), {'epsc', 'jpg', 'pdf'});
            % '_LD_causal_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
            
            
            figure; hold on; bar(non_causal_var_bins, non_causal_var_hist./sum(non_causal_var_hist), 0.2);
            bar(pruned_non_causal_var_bins+all_var_bin_size/4, pruned_non_causal_var_hist./sum(pruned_non_causal_var_hist), 0.2, 'r');
            bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist./sum(only_causal_var_hist), 0.2, 'g');
            xlabel('MAF'); ylabel('Var. Explained');
            legend({['All non-LD-causal (V=' num2str(V_non_causal,4) ')'], ...
                ['Pruned non-LD-causal (V=' num2str(V_pruned_non_causal,4) ')'], ...
                ['True (V=' num2str(V_true_FDR, 4) ' of ' num2str(V_true, 4) ')']}, 2);
            my_title(['Estimated Variance Explained, ' model_str]);
            xlim([0 0.51]); ylim(all_y_lim);
            my_saveas(gcf, fullfile(SNPs_figs_dir, [remove_dir_from_file_name(GWAS_results_file) 'non_LD_causal_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
        end % if plot by LD
        
        %            causal_LD_SNPs = cell2vec(R_LD.LD_SNPs(I_causal));
        
        %    hold on; % Plot true vs. observed effect sizes for causal SNPs
        subplot(2,2,3); % figure;   % New: Put all in the same file   % figure;
        hold on;
        plot(R_causal.beta(I_causal2), GWAS_RES.beta(J_causal2), 'r*');
        rand_snp_inds = randperm(num_snps); rand_snp_inds = rand_snp_inds(1:1000);
        plot(zeros(1, 1000), GWAS_RES.beta(rand_snp_inds), '.'); % plot a set of random non-causal snps
        plot( linspace(min(0, min(R_causal.beta)), max(R_causal.beta), 100), ...
            linspace(min(0, min(R_causal.beta)), max(R_causal.beta), 100), 'k'); % y=x
        xlabel('true-effect-size'); ylabel('observed-effect-size');
        legend({'causal', 'non-causal', 'y=x'}, 4, 'fontsize', 7);
        title(['True vs. observed \beta. Mean = ' num2str(mean(GWAS_RES.beta(J_causal2)), 3) ...
            ', std='  num2str(std(GWAS_RES.beta(J_causal2)), 3)]);
        %         my_saveas(gcf, fullfile(SNPs_figs_dir, [remove_dir_from_file_name(GWAS_results_file) ...
        %             '_allele_freq_var_dist']), {'epsc', 'jpg', 'pdf'});
        %     my_saveas(gcf, fullfile(SNPs_figs_dir, ...
        %         [remove_dir_from_file_name(GWAS_results_file) '_true_and_observed_effect_sizes']), ...
        %         {'epsc', 'jpg', 'pdf'});
        
        
        %    figure;
        subplot(2,2,4); hold on; % another way of plotting causal SNPS
        J_non_causal = setdiff(1:length(GWAS_RES.beta), J_causal2);
        hist_density(GWAS_RES.beta, 1000);  y_lim = get(gca, 'ylim'); y_height = diff(y_lim);
        plot(GWAS_RES.beta(J_causal2), zeros(length(J_causal2), 1)-0.01*y_height, 'go');
        plot(R_causal.beta, zeros(num_causal_snps, 1)+0.01*y_height, 'r*');
        xlabel('\beta'); ylabel('freq.');
        legend({'all-snps', 'causal ($\beta_{true}$)', 'causal($\hat{\beta}$)'}, 'fontsize', 7, 'interpreter', 'latex');
        title(['Distribution of effect sizes \beta. St.d.(\beta) = ' ...
            num2str(std(GWAS_RES.beta(J_non_causal)), 3)]);
        my_saveas(gcf, fullfile(SNPs_figs_dir, run_str, [remove_dir_from_file_name(GWAS_results_file) ...
            '_allele_freq_var_dist']), {'pdf'}); % save only once in order to not polute output directory
        %     my_saveas(gcf, fullfile(SNPs_figs_dir, ...
        %         [remove_dir_from_file_name(GWAS_results_file) '_effect_sizes_dist']), ...
        %         {'epsc', 'jpg', 'pdf'});
        
    end % if isdeployd
    
    
end % if exist causal SNPs file



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal function: preprocess_LD
% Save LD information in new format: snp and all it's ld snps in the same line
% Input:
% LD_file - file name for LD information (in input format: snp-a snp-b r2 ...)
% save_txt_flag - flag saying if to save results also in .txt format (in addition to .mat format)
%
function preprocess_LD(LD_file, save_txt_flag)

if(~exist(file_name_to_mat(LD_file), 'file'))
    R_LD = ReadDataFile(LD_file, []);
    save(file_name_to_mat(LD_file), 'R_LD');
else
    R_LD = load(file_name_to_mat(LD_file));
    if(isfield(R_LD, 'R_LD'))
        R_LD = R_LD.R_LD;
    end
end
tmp_SNP_A = R_LD.SNP_A; % duplicate pairs data
R_LD.SNP_A = [R_LD.SNP_A' R_LD.SNP_B']'; % duplicate to make sure no pair is lost
R_LD.SNP_B = [R_LD.SNP_B' tmp_SNP_A']'; % duplicate to make sure no pair is lost
R_LD.R2 = [R_LD.R2' R_LD.R2']';

%LD_SNPs = unique(R_LD.SNP_A); % union(unique(R_LD.SNP_A), unique(R_LD.SNP_B)); % list all SNPs

%[intersect_SNPs2 I2 J2] = intersect_all(R_SNPs.rsid, R_LD.SNP_A)

[unique_SNPs_A unique_inds_cell unique_num] = unique_with_inds(R_LD.SNP_A); % keep indices of all SNPs
[~, unique_I_A unique_J_A] = unique(R_LD.SNP_A); % keep indices of all SNPs

[unique_SNPs_B, unique_I_B, unique_J_B] = unique(R_LD.SNP_B); % keep indices of secon SNPs

num_snps = length(unique_SNPs_A); % nunber of unique SNPs
R_LD.rsid = unique_SNPs_A; % erase previous rsid if exist
R_LD.LD_inds = cell(num_snps,1);
R_LD.LD_SNPs = cell(num_snps,1);
R_LD.LD_R2 = cell(num_snps,1);


for i=1:num_snps % loop on all SNPs and get their 'neighbours'
    if(mod(i,1000) == 0)
        run_snp_i = i
    end
    R_LD.LD_inds{i} = unique_J_B(unique(unique_inds_cell{i})); % indices of SNP in LD with i-th SNP
    R_LD.LD_SNPs{i} = cell2vec(R_LD.rsid(R_LD.LD_inds{i}), ', ');
    R_LD.LD_R2{i} = R_LD.R2(R_LD.LD_inds{i});
end
R_LD = rmfield(R_LD, 'SNP_A'); R_LD = rmfield(R_LD, 'SNP_B');
save_R2 = R_LD.R2; R_LD.R2 = R_LD.LD_R2; R_LD = rmfield(R_LD, 'LD_R2');
R_LD.CHR_A = [R_LD.CHR_A' R_LD.CHR_B']';
R_LD.chr = R_LD.CHR_A(unique_I_A);
R_LD.BP_A = [R_LD.BP_A' R_LD.BP_B']';
R_LD.pos = R_LD.BP_A(unique_I_A);
R_LD = rmfield(R_LD, 'CHR_A'); R_LD = rmfield(R_LD, 'CHR_B');
R_LD = rmfield(R_LD, 'BP_A'); R_LD = rmfield(R_LD, 'BP_B');

[pos_sorted pos_sort_perm] = sort(R_LD.pos); % sort positions
R_LD = sort_struct_by_field(R_LD, 'pos'); % sort according to genomic position
inv_pos_sort_perm = inv_perm(pos_sort_perm); % change indices of SNPs in LD
for i=1:num_snps
    R_LD.LD_inds{i} = inv_pos_sort_perm(R_LD.LD_inds{i});
end
R_LD.sparse_mat = sparse(unique_J_A, unique_J_B, save_R2);     % New: save also all data in a sparse matrix form
R_LD.sparse_mat = R_LD.sparse_mat(pos_sort_perm,pos_sort_perm);  % save sparse form also sorted by position !!!

R_LD.block_vec = divide_matrix_to_blocks(R_LD.sparse_mat);

R_LD.inds = (1:num_snps)';
R_LD = orderfields(R_LD, {'inds', 'rsid', 'chr', 'pos', 'LD_inds', 'LD_SNPs', 'R2', 'sparse_mat', 'block_vec'});
save([remove_suffix_from_file_name(LD_file) '_preprocessed.mat'], 'R_LD'); % save in binary .mat format

if(save_txt_flag) % save as .txt
    for i=1:num_snps % write inds and LD nicely
        if(mod(i,1000) == 0)
            write_nicely_i = i
        end
        R_LD.LD_inds{i} = cell2vec(num2str_cell(num2cell( R_LD.LD_inds{i})), ', ');
        R_LD.R2{i} = cell2vec(num2str_cell(num2cell( R_LD.R2{i}), 4), ', ');
    end
    R_LD = rmfield(R_LD, {'sparse_mat', 'block_vec'}); % don't write these to .txt
    WriteDataFile(R_LD, [remove_suffix_from_file_name(LD_file) '_preprocessed.txt']); % save in .txt format
end


% % %         % No need to intersect lists
% % %         [intersect_SNPs I J] = intersect(R_SNPs.rsid, unique_SNPs_A); % intersect two lists
% % %
% % %         num_snps = length(R_SNPs.rsid);
% % %         R_SNPs.LD_inds = cell(num_snps,1); % indices of snps in LD with each SNP
% % %         R_SNPs.LD_SNPs = cell(num_snps,1); % which snps are in LD with each SNP
% % %         R_SNPs.LD_R2 = cell(num_snps,1); % level of LD for each SNP
% % %         for i=1:length(I) % loop on all SNPs and get their 'neighbours'
% % %             if(mod(i,1000) == 0)
% % %                 run_snp_i = i
% % %             end
% % %             %    cur_LD_ind = J(i);
% % %             [unique_SNPs_in_LD, unique_SNPs_in_LD_inds] = unique(R_LD.SNP_B(unique_inds_cell{J(i)}));
% % %             cur_R_LD = R_LD.R2(unique_inds_cell{J(i)});
% % %             cur_inds = I(unique_J(unique_inds_cell{J(i)}));
% % %             R_SNPs.LD_inds{I(i)} = cell2vec(num2str_cell(cur_inds(unique_SNPs_in_LD_inds)), ', '); % indices of snps in LD with each SNP
% % %             R_SNPs.LD_SNPs{I(i)} = cell2vec(unique_SNPs_in_LD, ', ');
% % %             R_SNPs.LD_R2{I(i)} = cell2vec(num2str_cell(cur_R_LD(unique_SNPs_in_LD_inds), 3), ', ');
% % %         end
% % %
% % %         R_SNPs.mark_vec = zeros(num_snps,1); % indicator vector: 1 - keep it, 0 - remove (due to LD)
% % %
% % %         P = 1:length(fields(R_SNPs)); P(end-3) = P(end); P(end-2:end) = [P(end)-3, P(end)-2, P(end)-1];
% % %         R_SNPs = orderfields(R_SNPs, P);
% % %         WriteDataFile(R_SNPs, [remove_suffix_from_file_name(SNPs_file) '_with_LD.txt']); % save in .txt format
% % %         save([remove_suffix_from_file_name(SNPs_file) '_with_LD.mat'], 'R_SNPs'); % save in binary .mat format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal function: prepare LD matrix (not working yet)
%
% Input:
% LD_file - file with LD information in the current form (prefix)
%
function prepare_LD_matrix(LD_file)

% % % % if(exist(file_name_to_mat(LD_file), 'file'))  % load R_LD
% % % %     R_LD = load(file_name_to_mat(LD_file), 'rsid', 'LD_inds', 'chr', 'R2', 'pos'); % , 'R2'); % no need for R2
% % % %     if(isempty(fields(R_LD)))
% % % %         R_LD = load(file_name_to_mat(LD_file));
% % % %     end
% % % %     if(isfield(R_LD, 'R_LD'))
% % % %         R_LD = R_LD.R_LD;
% % % %     end
% % % % else
% % % %     R_LD = ReadDataFile(LD_file);
% % % % end
for chr = 1:22 % load a separate LD file for each chromosome
    R_LD_chr = load(file_name_to_mat([ LD_file '_chr' num2str(chr)])); % load current chromosome
    if(isfield(R_LD_chr, 'R_LD'))
        R_LD_chr = R_LD_chr.R_LD;
    end
    chr_inds = find(R_LD.chr == chr); chr_num_snps = length(chr_inds); start_ind = min(chr_inds)-1;
    R_LD.sparse_mat{chr} = sparse(chr_num_snps,chr_num_snps);
    for j=1:chr_num_snps
        if(mod(j,1000)==0)
            j_is = j
        end
        R_LD.sparse_mat{chr}(j, R_LD.LD_inds{j} - start_ind) = 1;
    end
    [sorted_pos sort_perm] = sort(R_LD_chr.pos); % sort by position
    R_LD.sparse_mat{chr} = R_LD.sparse_mat{chr}(sort_perm,sort_perm);
end

save([LD_file '_genomic.mat'], 'R_LD');  % save united LD file




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal function: perform deconvolution
%
% Input:
% beta_vec - vector of observed effect sizes
% R_LD -
%
% Output:
% beta_deconvolved_vec - original effect sizes of SNPs
%
function beta_deconvolved_vec = LD_deconvolution(beta_vec, f_vec, R_LD) % perform deconvulution

C = posdefrnd(num_snps); % C = eye(num_snps); % correlation matrix
SNP_CORR_mat = corr_mat_to_IBD_mat(C, f_vec); % This is the r^2 matrix

SNP_beta_mat = (R_LD.R2_mat .* ( sqrt(f_vec .* (1-f_vec))' * (1./sqrt(f_vec .* (1-f_vec))) ))';
%    figure; plot(SNP_CORR_mat(:), SNP_beta_mat(:), '.')

beta_deconvolved_vec = beta_vec' / SNP_beta_mat'; % sqrtm(SNP_CORR_mat);



