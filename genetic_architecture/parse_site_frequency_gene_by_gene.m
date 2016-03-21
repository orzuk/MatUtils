% Fit selection and frequency for each gene using exome data 
%
% Input:
% spectrum_data_dir - input directory with frequency specturm data
% spectrum_data_file - file with sequencing data (changed to cell array for multiple populations!!)
% output_data_dir - where to save output files and figures
% GeneStruct - file/structure with genomic information on each gene
% MutationRateTable - total mutation rate per-gene for different mutation types
% MutationTypes - strings of each mutation type
% gene_prefix - run only on genes matching this prefix (used to divide labor)
%
function parse_site_frequency_gene_by_gene(spectrum_data_dir, spectrum_data_file, output_data_dir, ...
    GeneStruct, MutationRateTable, MutationTypes, gene_prefix) % Estimate s and alpha for each gene in the genome - how???

Assign24MammalsGlobalConstants;
if(ischar(GeneStruct))
    gene_struct_input_file = GeneStruct;
    if(~exist([remove_suffix_from_file_name(gene_struct_input_file) '_unique.mat'], 'file')); % Save unique
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

if(~exist([remove_suffix_from_file_name(gene_struct_input_file) '_unique.mat'], 'file')); % Save unique genes
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
        load_fields = [load_fields {'REF', 'ALT', 'aminoAcidChange', 'gene_by_allele_type_pos_list', 'gene_by_allele_type_inds_list'}];
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
    
    SiteFreqSpecStruct{k}.good_allele_inds = sort(SiteFreqSpecStruct{k}.good_allele_inds); % sort to fix mutation types ordering
    SiteFreqSpecStruct{k}.unique_genes = upper(SiteFreqSpecStruct{k}.unique_genes);
    SiteFreqSpecStruct{k}.allele_types_ind = zeros(size(SiteFreqSpecStruct{k}.allele_types));
    [~, I, J] = my_intersect(genome_types, vec2row(lower(strrep_cell(SiteFreqSpecStruct{k}.allele_types, '-', '_'))));
    setdiff( lower(strrep_cell(SiteFreqSpecStruct{k}.allele_types, '-', '_')), lower(empty_cell_to_empty_str(genome_types)) )
    for i=1:length(genome_types_synonyms)
        [~, tmp_I, tmp_J] = intersect(genome_types_synonyms{i}, lower(strrep_cell(SiteFreqSpecStruct{k}.allele_types, '-', '_')));
        if(~isempty(tmp_I))
            J = [J tmp_J];
            I = [I repmat(i, length(tmp_J), 1)];
        end
    end
    [J, UJ] = unique(J); I = I(UJ); % get rid of multiplicities
    SiteFreqSpecStruct{k}.allele_types_ind(J) = I;
    
    allele_types_ind = SiteFreqSpecStruct{k}.allele_types_ind;     
    save(fullfile(spectrum_data_dir, spectrum_data_file{k}), '-append', 'allele_types_ind'); % Save new values into file 
    
end % loop on populations

%for c = 'A':'Z'  % Make directories
%    my_mkdir(gene_dir);
%end
ExonsGeneStruct = load(gene_struct_input_file, 'chr_vec', 'pos_start_vec', 'pos_end_vec', 'strand', 'gene_names', 'sort_perm'); % load information on genes. 

for i=1:num_genes % loop on genes
    run_gene = i
    if(strmatch(upper(gene_prefix), upper(GeneStruct.gene_names{i})))
        gene_header = upper(GeneStruct.gene_names{i}) % (1:2));  % Print gene name 
        gene_dir = fullfile(output_data_dir, gene_header(1), GeneStruct.gene_names{i}); % here save all information on gene
        my_mkdir(gene_dir);
        internal_plot_gene_stats(gene_header, i, GeneStruct, ExonsGeneStruct, ...
            SiteFreqSpecStruct, MutationRateTable, MutationTypes, gene_dir);
    end
    if(mod(i, 10) == 0) % avoid having too many figures open
        close all;
    end
end






