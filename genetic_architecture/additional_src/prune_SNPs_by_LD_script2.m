% Prune table of SNP pvals based on LD structure
%
% Input:
% SNPs_file - file with data on all SNPs
% LD_file - file with pairwise LD data
% GWAS_results_file - file with p-values, FDR etc. for each SNP
% causal_SNPs_file - New: add a possibility of the true causal SNPs. This is used to remove stuff in LD with them
% operation_flag - 1: parse LD files, 0 - prune
% save_txt_flag - 1: save data also in .txt tab-delimited format (default, slow). 0 - save only in .mat format
%
function prune_SNPs_by_LD_script(SNPs_file, LD_file, GWAS_results_file, causal_SNPs_file, ...
    operation_flag, save_txt_flag, trait_type)

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
        'operation_flag - what to do. prune_snps (default) or preprocess_ld\n\n' ...
        'Output pruned file (for the default prune_snps operation) will be in the file: [GWAS_results-without-suffix ''_with_LD_pruned.txt'']\n\n' ...
        'Example: ./prune_SNPs_by_LD_script_exe ''''  ' ...
        '/seq/orzuk/common_disease_model/data/fdr/PLINK/CEU_chr1.ld_preprocessed.txt ' ...
        '/seq/orzuk/common_disease_model/data/fdr/chr1_20_causal_snps_bin_40_50.txt ' ...
        '/seq/orzuk/common_disease_model/data/fdr/causal_loci_chr1.txt prune_snps\n\n' ...
        'For this example, output will be in: /seq/orzuk/common_disease_model/data/fdr/chr1_20_causal_snps_bin_40_50_with_LD_pruned.txt'])
    return;
end

switch machine
    case PC
        SNPs_data_dir = '../../common_disease_model/data/fdr';
    case UNIX
        SNPs_data_dir = '/seq/orzuk/common_disease_model/data/fdr';
end
SNPs_figs_dir = fullfile(SNPs_data_dir, 'figs');
run_str = 'all'; % 'all'; % 'small'; % 'all';
% LD_thresh = 0.5;

if(~exist('operation_flag', 'var') || isempty(operation_flag)) % set default: compare LD
    operation_flag = 0;
end
if(~exist('save_txt_flag', 'var') || isempty(save_txt_flag)) % set default: compare LD
    save_txt_flag = 1;
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

switch operation_flag
    case {1, 'preprocess_ld'} % parse LD files. Do not need the SNP files
        %         if(~exist(file_name_to_mat(SNPs_file), 'file'))
        %             R_SNPs = ReadDataFile(SNPs_file);
        %             save(file_name_to_mat(SNPs_file), 'R_SNPs');
        %         else
        %             load(file_name_to_mat(SNPs_file));
        %         end
        sprintf('preprocess LD')
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
        
        [unique_SNPs_B, unique_I_B, unique_J_B] = unique(R_LD.SNP_B);
        
        num_snps = length(unique_SNPs_A); % nunber of unique SNPs
        R_LD.rsid = unique_SNPs_A; % erase previous rsid if exist
        R_LD.LD_inds = cell(num_snps,1);
        R_LD.LD_SNPs = cell(num_snps,1);
        R_LD.LD_R2 = cell(num_snps,1);
        
        for i=1:num_snps % loop on all SNPs and get their 'neighbours'
            if(mod(i,1000) == 0)
                run_snp_i = i
            end
            R_LD.LD_inds{i} = unique_J_B(unique(unique_inds_cell{i}));
            R_LD.LD_SNPs{i} = cell2vec(R_LD.rsid(R_LD.LD_inds{i}), ', ');
            R_LD.LD_R2{i} = R_LD.R2(R_LD.LD_inds{i});
        end
        R_LD = rmfield(R_LD, 'SNP_A'); R_LD = rmfield(R_LD, 'SNP_B');
        R_LD = rmfield(R_LD, 'R2'); R_LD.R2 = R_LD.LD_R2; R_LD = rmfield(R_LD, 'LD_R2');
        R_LD.CHR_A = [R_LD.CHR_A' R_LD.CHR_B']';
        R_LD.chr = R_LD.CHR_A(unique_I_A);
        R_LD.BP_A = [R_LD.BP_A' R_LD.BP_B']';
        R_LD.pos = R_LD.BP_A(unique_I_A);
        R_LD = rmfield(R_LD, 'CHR_A'); R_LD = rmfield(R_LD, 'CHR_B');
        R_LD = rmfield(R_LD, 'BP_A'); R_LD = rmfield(R_LD, 'BP_B');
        
        R_LD = sort_struct_by_field(R_LD, 'pos'); % sort according to genomic position
        R_LD.inds = (1:num_snps)';
        R_LD = orderfields(R_LD, {'inds', 'rsid', 'chr', 'pos', 'LD_inds', 'LD_SNPs', 'R2'});
        save([remove_suffix_from_file_name(LD_file) '_preprocessed.mat'], 'R_LD'); % save in binary .mat format
        
        if(save_txt_flag) % save as .txt
            for i=1:num_snps % write inds and LD nicely
                R_LD.LD_inds{i} = cell2vec(num2str_cell(num2cell( R_LD.LD_inds{i})), ', ');
                R_LD.R2{i} = cell2vec(num2str_cell(num2cell( R_LD.R2{i}), 4), ', ');
            end
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
        
        
    case {0, 'prune_snps'} % prune LD
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
            reading_output_file = file_name_to_mat(GWAS_results_file)
            GWAS_RES = ReadDataFile(GWAS_results_file, file_name_to_mat(GWAS_results_file));
            save(file_name_to_mat(GWAS_results_file), '-struct', 'GWAS_RES');
        end
        if(~isfield(GWAS_RES, 'FDR')) % set field name
            GWAS_RES.FDR = GWAS_RES.fdr;
        end
        GWAS_RES = my_rmfield(GWAS_RES, 'fdr');
        if(~isfield(GWAS_RES, 'pval')) % set field name
            GWAS_RES.pval = GWAS_RES.p;
        end
        GWAS_RES = my_rmfield(GWAS_RES, 'p');
        GWAS_RES = sort_struct_by_field(GWAS_RES, 'FDR', 'descend');
        num_snps = length(GWAS_RES.rsid);
        [intersect_SNPs I J] = intersect(R_LD.rsid, GWAS_RES.rsid);
        
        if(~isfield(GWAS_RES, 'f_vec')) % fill allele frequencies
            GWAS_RES.f_vec = zeros(length(GWAS_RES.rsid), 1);
            if(isfield(R_LD, 'f_vec'))
                GWAS_RES.f_vec(J) = R_LD.f_vec(I);
            else % need to load also SNPs file
                %            now_SNPs_file_is = SNPs_file
                if(~exist(file_name_to_mat(SNPs_file), 'file'))  % load R_SNPs
                    R_SNPs = ReadDataFile(SNPs_file, []);
                    save(file_name_to_mat(SNPs_file), 'R_LD');
                else
                    load(file_name_to_mat(SNPs_file));
                end
                [~, I_SNPs J_SNPs] = intersect(R_SNPs.rsid, GWAS_RES.rsid);
                GWAS_RES.f_vec(J_SNPs) = R_SNPs.all_maf(I_SNPs);
            end
        end % if isfield f_vec
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
            figure; hold on; bar(all_var_bins, all_var_hist);
            bar(pruned_var_bins, pruned_var_hist, 'r');
            xlabel('MAF'); ylabel('Var. Explained');
            legend({['All (V=' num2str(V_all, 4) ')'], ['Pruned (V=' num2str(V_pruned, 4) ')']}, 2);
            my_title(['Estimated Variance Explained, ' model_str]);
            all_y_lim = get(gca, 'ylim')
            my_saveas(gcf, fullfile(SNPs_figs_dir, [remove_dir_from_file_name(GWAS_results_file) '_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
        end
        if(exist(causal_SNPs_file, 'file') && (~isempty(causal_SNPs_file))) % seperate into stuff in LD and not from causal SNPs
            R_causal = ReadDataFile(causal_SNPs_file); % load causal snps
            num_causal_snps = length(R_causal.rsid)
            [causal_SNPs I_causal J_causal] = intersect(R_LD.rsid, R_causal.rsid); % get indices of causal SNPs
            [~, I_causal2, J_causal2] = intersect(R_causal.rsid, GWAS_RES.rsid);
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
            if ~isdeployed
                bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist, 'g');
                legend({['All (V=' num2str(V_all, 4) ')'], ...
                    ['Pruned (V=' num2str(V_pruned, 4) ')'], ...
                    ['True (V=' num2str(V_true_FDR, 4) ' of ' num2str(V_true, 4) ')']}, 2);
                xlim([0 0.51]);
                my_saveas(gcf, fullfile(SNPs_figs_dir, ...
                    [remove_dir_from_file_name(GWAS_results_file) '_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
            end
            
            [~, I_causal2, causal_LD_SNP_inds2] = intersect(R_LD.rsid(causal_LD_SNP_inds), GWAS_RES.rsid);
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
                figure; hold on; bar(all_var_bins, all_var_hist);
                bar(pruned_var_bins, pruned_var_hist, 'r');
                xlabel('MAF'); ylabel('Var. Explained');
                legend({['All (V=' num2str(V_all, 4) ')'], ['Pruned (V=' num2str(V_pruned, 4) ')']}, 2);
                my_title(['Estimated Variance Explained, ' model_str ', p-val < ' num2str(pval_cutoff, 3)]);
                bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist, 'g');
                legend({['All (V=' num2str(V_all, 4) ')'], ...
                    ['Pruned (V=' num2str(V_pruned, 4) ')'], ...
                    ['True (V=' num2str(V_true_FDR, 4) ' of ' num2str(V_true, 4) ')']}, 2);
                xlim([0 0.51]); ylim(all_y_lim);
                my_saveas(gcf, fullfile(SNPs_figs_dir, ...
                    [remove_dir_from_file_name(GWAS_results_file) '_var_explained_with_pruning_high_res_strict_pval']), {'epsc', 'jpg', 'pdf'});
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
                figure; hold on; bar(causal_var_bins, causal_var_hist);
                bar(pruned_causal_var_bins, pruned_causal_var_hist, 'r');
                bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist, 'g');
                xlabel('MAF'); ylabel('Var. Explained');
                legend({['All LD-causal (V=' num2str(V_causal,4) ')'], ...
                    ['Pruned LD-causal (V=' num2str(V_pruned_causal,4) ')'], ...
                    ['True (V=' num2str(V_true_FDR, 4) ' of ' num2str(V_true, 4) ')']}, 2);
                my_title(['Estimated Variance Explained, ' model_str]);
                xlim([0 0.51]); ylim(all_y_lim);
                my_saveas(gcf, fullfile(SNPs_figs_dir, [remove_dir_from_file_name(GWAS_results_file) '_LD_causal_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
                
                
                figure; hold on; bar(non_causal_var_bins, non_causal_var_hist);
                bar(pruned_non_causal_var_bins, pruned_non_causal_var_hist, 'r');
                bar(only_causal_var_bins+only_causal_bin_size/2, only_causal_var_hist, 'g');
                xlabel('MAF'); ylabel('Var. Explained');
                legend({['All non-LD-causal (V=' num2str(V_non_causal,4) ')'], ...
                    ['Pruned non-LD-causal (V=' num2str(V_pruned_non_causal,4) ')'], ...
                    ['True (V=' num2str(V_true_FDR, 4) ' of ' num2str(V_true, 4) ')']}, 2);
                my_title(['Estimated Variance Explained, ' model_str]);
                xlim([0 0.51]); ylim(all_y_lim);
                my_saveas(gcf, fullfile(SNPs_figs_dir, [remove_dir_from_file_name(GWAS_results_file) 'non_LD_causal_var_explained_with_pruning_high_res']), {'epsc', 'jpg', 'pdf'});
            end
            
            
            %            causal_LD_SNPs = cell2vec(R_LD.LD_SNPs(I_causal));
            
            
        end
        
    case 2 % prepare deconvolution matrices 
        if(exist(file_name_to_mat(LD_file), 'file'))  % load R_LD
            R_LD = load(file_name_to_mat(LD_file), 'rsid', 'LD_inds', 'chr', 'R2', 'pos'); % , 'R2'); % no need for R2
            if(isempty(fields(R_LD)))
                R_LD = load(file_name_to_mat(LD_file));
            end
            if(isfield(R_LD, 'R_LD'))
                R_LD = R_LD.R_LD;
            end
        else
            R_LD = ReadDataFile(LD_file);
        end
        for chr = 1:22
            R_LD_chr = load(file_name_to_mat(LD_file));
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
               sort_perm = sort(R_LD_chr.pos); % sort by position 
           end
        end
        
end % switch operation_flag



% Internal function: Plot var explained

