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
%               There is also a file for each gene saved at: spectrum_data_dir
%               which contains fitted parameters 
%
function parse_site_frequency_gene_by_gene(spectrum_data_dir, spectrum_data_file, output_data_dir, ...
    GeneStruct, MutationRateTable, MutationTypes, Demographic_model, plot_flag, gene_prefix) % Estimate s and alpha for each gene in the genome - how???

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
if(ischar(MutationRateTable))
    MutationRateTable_file = MutationRateTable;
    load(MutationRateTable); % , 'MutationRateTable', 'MutationTypes');    
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
    save(MutationRateTable_file, '-append', 'UniqueMutationRateTable'); % save unique 
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

num_populations = length(spectrum_data_file);
SiteFreqSpecStruct = cell(num_populations, 1);
for k=1:num_populations % load data from all populations
    load_fields = {'unique_genes', 'n_vec', 'count_vec', 'f_vec', 'allele_types', ...
        'num_allele_types', 'num_alleles_per_gene_mat', 'total_freq_per_gene_mat', 'total_heterozygosity_per_gene_mat', ...
        'gene_by_allele_type_freq_list', 'gene_by_allele_type_n_list', ...
        'gene_by_allele_type_het_list', ...
        'good_allele_inds', 'upper_freq_vec'};
    if(k == 1)  % first population 
        load_fields = [load_fields {'REF', 'ALT', 'aminoAcidChange', 'gene_by_allele_type_pos_list', 'gene_by_allele_type_inds_list', 'allele_types_ind'}];
    end
    load_fields_str = cell2vec(load_fields, ''', ''');
    load_str = ['SiteFreqSpecStruct{' num2str(k) '} = load(''' fullfile(spectrum_data_dir, spectrum_data_file{k}) ...
        ''', ''' load_fields_str ''');'];
    eval(load_str);
%    SiteFreqSpecStruct{k} = load(fullfile(spectrum_data_dir, spectrum_data_file{k}), load_fields);
    % % % %         'unique_genes', 'n_vec', 'count_vec', 'f_vec', 'allele_types', ...
    % % % %         'num_allele_types', 'num_alleles_per_gene_mat', 'total_freq_per_gene_mat', 'total_heterozygosity_per_gene_mat', ...
    % % % %         'gene_by_allele_type_freq_list', 'gene_by_allele_type_n_list', ...
    % % % %         'gene_by_allele_type_het_list', 'gene_by_allele_type_pos_list', 'gene_by_allele_type_inds_list', ...
    % % % %         'good_allele_inds', 'upper_freq_vec', ...
    % % % %         'REF', 'ALT', 'aminoAcidChange'); %
    %        'het_vec', 'het_var_vec', 'variants', 'carriers', 'singletons', 'heterozygosity'); % load sequencing input data from file
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
    end
        
end % loop on populations

ExonsGeneStruct = load(gene_struct_input_file, ...
    'chr_vec', 'pos_start_vec', 'pos_end_vec', 'strand', 'gene_names', 'sort_perm'); % load information on genes. 

for i=1:num_genes % loop on genes and plot / fit selection coefficients 
    sprintf(['Run gene = %d out of %d, ' upper(GeneStruct.gene_names{i})], i, num_genes)
    if(startsWith(upper(GeneStruct.gene_names{i}), upper(gene_prefix)))
        gene_header = upper(GeneStruct.gene_names{i}); % (1:2));  % Print gene name 
        gene_dir = fullfile(output_data_dir, gene_header(1), GeneStruct.gene_names{i}); % here save all information on gene
        my_mkdir(gene_dir);
        if(plot_flag)
            internal_plot_gene_stats(gene_header, i, GeneStruct, ExonsGeneStruct, ...
                SiteFreqSpecStruct, MutationRateTable, MutationTypes, gene_dir, plot_flag);
        end
        if(fit_selection)
            alpha_vec = 0.1:0.1:1; % possible values for alpha 
            s_null_vec = 0 -logspace(-6, -1, 10);
            rare_cumulative_per_gene = 1; 
            target_size_by_class_vec = [mu_per_site, mu_per_site, mu_per_site]; % Get target size for gene
            gene_n_vec = []; gene_k_vec = []; null_w_vec = []; % ???
            for j= SiteFreqSpecStruct{1}.good_allele_inds{2} % loop on missense
                gene_k_vec = [gene_k_vec (SiteFreqSpecStruct{1}.gene_by_allele_type_freq_list{j,i} .* ...
                    SiteFreqSpecStruct{1}.gene_by_allele_type_n_list{j,i})'];
                gene_n_vec = [gene_n_vec SiteFreqSpecStruct{1}.gene_by_allele_type_n_list{j,i}']; 
            end
            null_w_vec = -ones(size(gene_k_vec)); % all missense
            for j= SiteFreqSpecStruct{1}.good_allele_inds{3} % loop on stop codons 
                gene_k_vec = [gene_k_vec (SiteFreqSpecStruct{1}.gene_by_allele_type_freq_list{j,i} .* ...
                    SiteFreqSpecStruct{1}.gene_by_allele_type_n_list{j,i})'];
                gene_n_vec = [gene_n_vec SiteFreqSpecStruct{1}.gene_by_allele_type_n_list{j,i}']; 
            end
            X = [gene_k_vec gene_n_vec]';     
            null_w_vec((end+1):length(gene_k_vec)) = 1; % all stop

            num_individuals = [];
            implementation_str = []; 
            trait_struct = []; % doesn't matter, we only use genotype part 
            maximize_parameters = [1 1 0]; % maximize over s and alpa 
            full_flag = 0; 
            if(~isempty(X)) % for some genes we have no alleles in the relevant population
                for k=2 % 1:num_populations % load data from all populations
                    [max_LL_vec(i,k), s_MLE_vec(i,k), alpha_MLE_vec(i,k), ~, max_compute_time(i,k)] = ... % don't fit  beta_MLE_vec(i,k)] = ...
                        maximize_two_class_likelihood(s_null_vec, alpha_vec, 0, ...
                        target_size_by_class_vec, Demographic_model{k}, ... % run each time on different population. But counts should be different!!
                        X, [], trait_struct, null_w_vec, maximize_parameters, full_flag, ...
                        num_individuals, implementation_str); % fit for each population seperately
                end % loop on populations
            end
        end % if fit selection
        
        aggregate_population_estimators = 0; 
        if(aggregate_population_estimators) % compute an aggregate esitmator for each gene from multiple populations
        
        
        end
        test_population_differences = 0; 
        if(test_population_differences) % test for differences in selection coefficient between different populations
            for population1 = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
                for population2 = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
                    
                    
                end
            end
        end
            
    end % filter for genes 
    if(mod(i, 10) == 0) % avoid having too many figures open
        close all;
    end
end % loop on genes






