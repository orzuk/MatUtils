% Fit selection and frequency for each gene using exome data
%
% Input:
% spectrum_data_dir - input directory with frequency specturm data
% spectrum_data_file - file with sequencing data (changed to cell array for multiple populations!!)
% output_data_dir - where to save output files and figures
% GeneStruct - file/structure with genomic information on each gene
% MutationRateTable - total mutation rate per-gene for different mutation types
% MutationTypes - strings of each mutation type
% Demographic_model - model for demography for each popoulation
% plot_flag - 1 (default): make plots, 0 - don't plot just save files
% gene_prefix - run only on genes matching this prefix (used to divide labor, run parallel jobs)
%
% Output: None.
%               Figures are saved in directory: output_data_dir
%               There is also a summary file with fitted parameters at: spectrum_data_dir
%               There is also a file for each gene saved at: spectrum_data_dir which contains fitted parameters
%
function parse_site_frequency_gene_by_gene(spectrum_data_dir, spectrum_data_file, output_data_dir, ...
    GeneStruct, MutationRateTable, MutationTypes, Demographic_model, plot_flag, gene_prefix) % Estimate s and alpha for each gene in the genome - how???

if(ischar(spectrum_data_file))
    spectrum_data_file = strcat(dir_from_file_name(spectrum_data_file), ...
        GetFileNames(fullfile(spectrum_data_dir, spectrum_data_file))); % here for all files in chunks
end
num_files = length(spectrum_data_file);
exome_data = str2word('\', spectrum_data_file(1), 1); exome_struct = get_exome_data_info(exome_data{1});
Assign24MammalsGlobalConstants; AssignRVASConstants;
if(~exist('plot_flag', 'var') || isempty(plot_flag))
    plot_flag = 0;
end
if(~exist('fit_selection', 'var') || isempty(fit_selection))
    fit_selection = 1; % fit selection coefficient for each gene
end
if(ischar(GeneStruct))
    gene_struct_input_file = GeneStruct;
    if(~exist([remove_suffix_from_file_name(gene_struct_input_file) '_unique.mat'], 'file')) % Save unique
        GeneStruct = load(GeneStruct, ...
            'chr_vec', 'pos_start_vec', 'pos_end_vec', 'seqs', 'strand', 'gene_names', 'sort_perm');
    end
end
if(ischar(MutationRateTable)) % load mutation table
    MutationRateTable_file = MutationRateTable;
    MutStruct = load(MutationRateTable); % , 'MutationRateTable', 'MutationTypes');
    if(isfield(MutStruct, 'MutationRateTable'))
        MutationRateTable = MutStruct.MutationRateTable;
    else % parse Samocha et al.'s file
        MutationRateTable = zeros(length(MutStruct.gene), 4);
        MutationRateTable(:,1) = 10.^MutStruct.syn;
        MutationRateTable(:,2) = 10.^MutStruct.mis;
        MutationRateTable(:,3) = 10.^MutStruct.non;
        MutationRateTable(:,4) = 10.^MutStruct.frameshift;
        MutStruct.MutationTypes = [SYNONYMOUS MISSENSE STOP STOP_LOST];
        MutStruct.gene_names = MutStruct.gene; % copt genes
        MutStruct.MutationRateTable = MutationRateTable; % copy to structure
    end
end
if(~exist('gene_prefix', 'var'))
    gene_prefix = ''; % empty string to match everything
end

num_mutation_types = length(MutationTypes);
TotalExonicMutationRateVec = sum(MutationRateTable); % Compute total mutation rate for each mutation to get target size

if(~exist([remove_suffix_from_file_name(gene_struct_input_file) '_unique.mat'], 'file')) % Save unique genes
    [UniqueGeneStruct.gene_names, I, J] = unique(GeneStruct.gene_names);  % Unite exons into genes
    num_genes = length(I); % Get only unique genes
    UniqueMutationRateTable = zeros(num_genes, num_mutation_types);
    for i=1:num_mutation_types
        UniqueMutationRateTable(:,i) = accumarray(J, MutationRateTable(:,i));
    end
    gene_names = UniqueGeneStruct.gene_names;
    save(MutationRateTable_file, '-append', 'UniqueMutationRateTable', 'gene_names'); % save unique mutation rates, with gene names (matching unique GeneStruct)
    MutationRateTable = UniqueMutationRateTable; % copy table of unique ones
    UniqueGeneStruct.chr_vec = GeneStruct.chr_vec(I);
    UniqueGeneStruct.seqs = cell(num_genes, 1);   % Unite sequences
    for i=1:num_genes
        if(mod(i, 500) == 0)
            unite_gene = i
        end
        UniqueGeneStruct.seqs{i} = cell2vec(GeneStruct.seqs(J == i));
    end
    UniqueGeneStruct.gene_lens = length_cell(UniqueGeneStruct.seqs);
    UniqueGeneStruct.strand = GeneStruct.strand(I);
    GeneStruct = UniqueGeneStruct; GeneStruct.gene_names = upper(GeneStruct.gene_names); % copy unique
    save([remove_suffix_from_file_name(gene_struct_input_file) '_unique.mat'], 'GeneStruct'); % Save unique
else % load unique genes
    load([remove_suffix_from_file_name(gene_struct_input_file) '_unique.mat']); % , ...
    %        'chr_vec', 'pos_start_vec', 'pos_end_vec', 'seqs', 'strand', 'gene_names', 'sort_perm');
end
if(exist('UniqueMutationRateTable', 'var'))
    MutationRateTable = UniqueMutationRateTable;
end
num_genes = length(GeneStruct.gene_names); % get total # of genes
num_populations = length(exome_struct.pop_str);

SiteFreqSpecStruct = cell(num_files, 1);
%for k=1:num_populations % load data from all populations%
%end % loop on populations

ExonsGeneStruct = load(gene_struct_input_file, ...
    'chr_vec', 'pos_start_vec', 'pos_end_vec', 'strand', 'gene_names', 'sort_perm'); % load information on genes.
[GeneStruct.s_MLE_vec, GeneStruct.s_MLE_vec, GeneStruct.alpha_MLE_vec, GeneStruct.max_compute_time] = ...
    deal(zeros(num_genes, num_populations));

load_fields = {'unique_genes', 'n_vec', 'count_vec', 'f_vec', 'allele_types', ...
    'num_allele_types', 'num_alleles_per_gene_mat', 'total_freq_per_gene_mat', 'total_heterozygosity_per_gene_mat', ...
    'gene_by_allele_type_freq_list', 'gene_by_allele_type_n_list', ...
    'gene_by_allele_type_het_list', ...
    'good_allele_inds', 'upper_freq_vec'};
%    if(k == 1)  % first population
load_fields = [load_fields {'REF', 'ALT', 'aminoAcidChange', 'gene_by_allele_type_pos_list', 'gene_by_allele_type_inds_list', 'allele_types_ind'}];
%    end
all_fit_genes_I = [];
for k=1:num_files % loop on different chunks
    sprintf('Run genes file %d out of %d', k, num_files) % plot progress
    load_fields_str = cell2vec(load_fields, ''', ''');
    load_str = ['SiteFreqSpecStruct{' num2str(k) '} = load(''' fullfile(spectrum_data_dir, spectrum_data_file{k}) ...
        ''', ''' load_fields_str ''');'];
    eval(load_str);
    if(~isfield(SiteFreqSpecStruct{k}, 'population_str'))
        SiteFreqSpecStruct{k}.population_str = ...
            str2word('_', remove_suffix_from_file_name(remove_dir_from_file_name(spectrum_data_file{k})), 'end');
    end
    if(~isfield(SiteFreqSpecStruct{k}, 'gene_by_allele_type_het_list'))
        for i=1:num_mutation_types
            for j=1:num_genes
                SiteFreqSpecStruct{k}.gene_by_allele_type_het_list{i,j} = 2 .* ...
                    SiteFreqSpecStruct{k}.gene_by_allele_type_freq_list{i,j} .* ...
                    (1 - SiteFreqSpecStruct{k}.gene_by_allele_type_freq_list{i,j}); % get heterozygosity
            end
        end
        gene_by_allele_type_het_list = SiteFreqSpecStruct{k}.gene_by_allele_type_freq_list;
        save('-append', fullfile(spectrum_data_dir, spectrum_data_file{k}), 'gene_by_allele_type_het_list');
        %        clear gene_by_allele_type_freq_list;
    end   % finished loading SFS
    
    [fit_genes, fit_genes_I, fit_genes_J] = ...
        intersect( upper(GeneStruct.gene_names), upper(SiteFreqSpecStruct{k}.unique_genes));
    all_fit_genes_I = [all_fit_genes_I fit_genes_I'];
    ctr=1;
    
    if(fit_selection)
        % First set parameters:
        alpha_vec = 0.1:0.1:1; % possible values for alpha
        s_null_vec = [0 -logspace(-6, -1, 10)];
        beta_vec = 0;
        target_size_by_class_vec = [mu_per_site, mu_per_site, mu_per_site]; % Get target size for gene
        num_individuals = []; implementation_str = []; trait_struct = []; % doesn't matter, we only use genotype part
        maximize_parameters = [1 1 0]; % maximize over s and alpha
        fit_type = 'total_load'; % 'MLE'; % total_load
        full_flag = 0;
        [GeneStruct.max_LL_vec(fit_genes_I,:), GeneStruct.s_MLE_vec(fit_genes_I,:), GeneStruct.alpha_MLE_vec(fit_genes_I,:), ~, GeneStruct.max_compute_time(fit_genes_I,:)] = ...
            fit_selection_parameters_internal(GeneStruct, ExonsGeneStruct, SiteFreqSpecStruct{k}, MutStruct, fit_genes_I, fit_genes_J, ...
            trait_struct, maximize_parameters, full_flag, fit_type, ...
            s_null_vec, alpha_vec, 0, ...
            target_size_by_class_vec, Demographic_model, num_populations, plot_flag, gene_prefix, ... % run each time on different population. But counts should be different!!
            num_individuals, implementation_str, output_data_dir); %
    end
    
end % loop on chunks (different files)
internal_save_gene_stats(GeneStruct, ExonsGeneStruct, output_data_dir, gene_struct_input_file, exome_struct, all_fit_genes_I);   % save fitted values to .mat and .txt files


% loop over genes and fit selection coefficient
function [LL_tab, s_hat, alpha_hat, beta_hat, compute_time] = ...
    fit_selection_parameters_internal(GeneStruct, ExonsGeneStruct, SiteFreqSpecStruct, MutStruct, fit_genes_I, fit_genes_J, ...
    trait_struct, maximize_parameters, full_flag, fit_type, ...
    s_null_vec, alpha_vec, beta_vec, ...
    target_size_by_class_vec, Demographic_model, num_populations, plot_flag, gene_prefix, ... % run each time on different population. But counts should be different!!
    num_individuals, implementation_str, output_data_dir) %

Assign24MammalsGlobalConstants; AssignRVASConstants;
[LL_tab,  s_hat, alpha_hat, beta_hat, compute_time] = deal(zeros(length(fit_genes_I), num_populations));

i_ctr=1;
for i=vec2row(fit_genes_I) % 1:num_genes % loop on genes and plot / fit selection coefficients
    if(plot_flag)
        sprintf(['Run gene = %d out of %d, ' upper(GeneStruct.gene_names{i})], i_ctr, length(fit_genes_I))
    end
    if(startsWith(upper(GeneStruct.gene_names{i}), upper(gene_prefix)))
        gene_header = upper(GeneStruct.gene_names{i}); % (1:2));  % Print gene name
        gene_dir = fullfile(output_data_dir, gene_header(1), GeneStruct.gene_names{i}); % here save all information on gene
        my_mkdir(gene_dir);
        if(plot_flag)
            internal_plot_gene_stats(gene_header, i, GeneStruct, ExonsGeneStruct, ...
                SiteFreqSpecStruct, MutStruct.MutationRateTable, MutStruct.MutationTypes, gene_dir, plot_flag);
        end
        % find gene's SFS
        i_sfs = fit_genes_J(i_ctr); %ctr=ctr+1;
        rare_cumulative_per_gene = 1;
        gene_n_vec = []; gene_k_vec = []; null_w_vec = []; % ???
        for j= SiteFreqSpecStruct.good_allele_inds{2} % loop on missense
            gene_k_vec = [gene_k_vec (SiteFreqSpecStruct.gene_by_allele_type_freq_list{j,i_sfs} .* ...
                SiteFreqSpecStruct.gene_by_allele_type_n_list{j,i_sfs})'];
            gene_n_vec = [gene_n_vec SiteFreqSpecStruct.gene_by_allele_type_n_list{j,i_sfs}'];
        end
        null_w_vec = -ones(1,size(gene_k_vec,2)); % all missense
        for j= SiteFreqSpecStruct.good_allele_inds{3} % loop on stop codons
            gene_k_vec = [gene_k_vec (SiteFreqSpecStruct.gene_by_allele_type_freq_list{j,i_sfs} .* ...
                SiteFreqSpecStruct.gene_by_allele_type_n_list{j,i_sfs})'];
            gene_n_vec = [gene_n_vec SiteFreqSpecStruct.gene_by_allele_type_n_list{j,i_sfs}'];
        end
        null_w_vec((end+1):length(gene_k_vec)) = 1; % all stop
        
        
        if(~isempty(gene_n_vec)) % for some genes we have no alleles in the relevant population
            for i_pop=1:num_populations % load data from all populations
                gene_n_vec_pop = gene_n_vec(i_pop, ~isnan(gene_k_vec(i_pop,:))); % Remove nan cells
                gene_k_vec_pop = gene_k_vec(i_pop, ~isnan(gene_k_vec(i_pop,:)));
                null_w_vec_pop = null_w_vec(~isnan(gene_k_vec(i_pop,:)));
                X = round([gene_k_vec_pop gene_n_vec_pop]');
                if(isempty(gene_n_vec_pop) || (max(gene_n_vec_pop) == 0))
                    continue; % go to next population
                end
                switch fit_type
                    case 'MLE'
                        [LL_tab(i_ctr,i_pop), s_hat(i_ctr,i_pop), alpha_hat(i_ctr,i_pop), ~, compute_time(i_ctr,i_pop)] = ... % don't fit  beta_MLE_vec(i,i_pop)] = ...
                            maximize_two_class_likelihood(s_null_vec, alpha_vec, beta_vec, ...
                            target_size_by_class_vec, Demographic_model{i_pop}, ... % run each time on different population. But counts should be different!!
                            X, [], trait_struct, null_w_vec, maximize_parameters, full_flag, ...
                            num_individuals, implementation_str); % fit for each population seperately
                    case 'total_load'
                        i_mut = strmatch(GeneStruct.gene_names{i}, MutStruct.gene_names, 'exact');
                        if(~isempty(i_mut))
                            f_lof = sum(gene_k_vec_pop(null_w_vec_pop==1) ./ gene_n_vec_pop(null_w_vec_pop==1));
                            f_missense = sum(gene_k_vec_pop(null_w_vec_pop==-1) ./ gene_n_vec_pop(null_w_vec_pop==-1));
                            if(f_lof==0) % no LOF varaints at all! complete depletion
                                s_hat(i_ctr,i_pop) = 1;
                            else
                                s_hat(i_ctr,i_pop) = sum(MutStruct.MutationRateTable(i_mut,((MutStruct.MutationTypes == STOP) ...
                                    | (MutStruct.MutationTypes == STOP_LOST)))) / f_lof; % take s_hat = mu / f_LOF
                            end
                            if(isnan(s_hat(i_ctr,i_pop)))
                                nananana = 999999999999
                            end
                            alpha_hat(i_ctr,i_pop) = sum(MutStruct.MutationRateTable(i_mut,((MutStruct.MutationTypes == MISSENSE)))) ...
                                ./ f_missense; % take ratio of expected loads of missense vs. stop codons
                        end
                end
                if(plot_flag)
                    sprintf('s = %f, alpha=%f, LL=%f, run-time (sec.)=%f, ', ...
                        s_hat(i_ctr,i_pop), alpha_hat(i_ctr,i_pop), ...
                        LL_tab(i_ctr,i_pop), compute_time(i_ctr,i_pop))
                end
            end % loop on populations
        end
        % Below add combining and testing for different populations
        % % %         aggregate_population_estimators = 0;
        % % %         if(aggregate_population_estimators) % compute an aggregate esitmator for each gene from multiple populations
        % % %         end
        % % %         test_population_differences = 0;
        % % %         if(test_population_differences) % test for differences in selection coefficient between different populations
        % % %             for population1 = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
        % % %                 for population2 = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
        % % %                 end
        % % %             end
        % % %         end
    end % filter for genes
    if(mod(i, 10) == 0) % avoid having too many figures open
        close all;
    end
    i_ctr=i_ctr+1;
end % loop on genes in chunk (file)



