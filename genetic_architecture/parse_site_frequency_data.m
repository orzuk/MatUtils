% Read and process the data file containing site-frequency-spectrum
%
% Input:
% site_frequency_file_name - file with information on each allele
% gene_list - list of possible gene names
% read_to_mat_flag - initial parsing. Read tab-delimited file with ... and convert it to .mat file (to fill!)
% extract_fields_flag - second parsing. Fill in field with new field names ... this convert the text extracted from text file to numerical values in matlab vectors (to fill!)
% compute_gene_matrices_flag - third parsing. Collect set of variants for each gene
%
% Output:
% S - structure with information on each allele
% n_vec - # of individuals with genotypes in each class
% count_vec -  # of individuals with alternative allele in each class
% f_vec - allele frequency of alternative variants in each class
%
function [S, n_vec, count_vec, f_vec] = parse_site_frequency_data(site_frequency_file_name, gene_list, ...
    read_to_mat_flag, extract_fields_flag, compute_gene_matrices_flag)

S = []; 
Assign24MammalsGlobalConstants; % get some constants needed
populations_vec = {''}; % {'_European', '_African'};
compute_frac_carriers = 1; % compute how many carriers of alleles below a certain frequency
compute_selection = 0; % fit selection coefficient.

if(~compute_gene_matrices_flag) % NEW! don't return if computing matrices!!! This should be done only once on the united file!!
    if(  (~exist(site_frequency_file_name, 'file')) && (~exist(file_name_to_mat(site_frequency_file_name), 'file')) )
        S = []; return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Stage 1: Convert to .mat file.           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return_flag=0
if(read_to_mat_flag) % convert vcf to mat files
    if(~exist(file_name_to_mat(site_frequency_file_name), 'file')) % read input .vcf file and convert to .mat
        reading_vcf_file=1
        %    populations_vec = {'_European', '_African'};
        if(strfind(lower(site_frequency_file_name), 'vcf'))
            S = read_vcf(site_frequency_file_name, file_name_to_mat(site_frequency_file_name)); % read as vcf file
        else
            S = ReadDataFile(site_frequency_file_name, file_name_to_mat(site_frequency_file_name)); % , cell_to_mat, skip_lines, delimiter, varargin)
        end
    else % read already made .mat file
        load_file_mat_format=1
        S = load(file_name_to_mat(site_frequency_file_name));
        
        preprocess_mat=1
        if(isfield(S, 'POS') && iscell(S.POS)) % re-format specific fields to save memory and time
            S.POS = cell2mat(str2num_cell(S.POS));
        end
        if(isfield(S, 'XXX_CHROM') && iscell(S.XXX_CHROM)) % re-format specific fields to save memory and time
            S.XXX_CHROM = chr_str2num(S.XXX_CHROM);
        end
        
        S.field_descriptions = strrep_cell(S.field_descriptions, 'Esitmated', 'Estimated'); % re-format descriptions to make them unique !!!
        [dup_vals, dup_inds, dup_num] = get_duplicates(S.field_descriptions);
        for j=1:length(dup_num) % change field names for duplicates
            if(dup_num(j) > 1) % duplicates
                for k=dup_inds{j}
                    S.field_descriptions{k} = [S.field_descriptions{k} '_' S.field_names{k}];
                end
            end
        end
    end
    S = my_rmfield(S, 'QUAL'); S = my_rmfield(S, 'FILTER'); S = my_rmfield(S, 'CoverageMat'); % Get rid of redundant fields
    S.field_names = S.field_names(1,:);
    S.field_descriptions = S.field_descriptions(1,:);
    save_file_mat_format=1
    save(file_name_to_mat(site_frequency_file_name), '-struct', 'S'); % always save in new format
    
end % if read_to_mat_flag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% END Stage 1: Convert to .mat file.         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Stage 2: Extract All Info Fields to Population file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compute_snp_info=1
if( strcmp(suffix_from_file_name(site_frequency_file_name), 'mat') )
    %    populations_vec = {''};
    populations_vec = {'_European', '_African'}; % special for ESP populations
else
    populations_vec = {'_European', '_African'};
    return_flag = 1
end
if(extract_fields_flag)
    for population = populations_vec  % save different files per population
        if(exist([remove_suffix_from_file_name(site_frequency_file_name) population{1} '.mat'], 'file')) % take only valuse which appear
            continue;
        end
        
        % Further parsing of vcf file
        num_snps = length(S.INFO);
        num_fields = length(S.field_names);
        S.XXX_VARIANT_COUNT_= zeros(num_snps,1);
        S.XXX_REF_ALLELE_COUNT_= zeros(num_snps,1);
        S.XXX_FEATURE_ = cell(num_snps,1);
        S.GENE = cell(num_snps,1);
        S.INFO_ARR = cell(num_snps, num_fields); % create big array of all fields !!
        %%%%        for j=1:length(S.field_names) % loop on all fields - this may be unnecessary
        %%%%            prepare_field=j
        %%%%            eval_str = ['S.' S.field_descriptions{j} ' = cell(' num2str(num_snps) ', 1);'];
        %%%%            eval(eval_str);
        %%%%        end
        start_snp_loop_time = cputime;
        for i=1:num_snps  % heavy loop: run on all SNPs and parse information
            if(mod(i, 500) == 0)
                sprintf('Parse SNP %d out of %d, time=%f', i, num_snps, cputime - start_snp_loop_time)
            end
            tmp_str = regexp(S.INFO{i}, ';', 'split');
            switch population{1}
                case {'European', '_European'}
                    tmp_allele_counts = str2nums(tmp_str{2}); % Take Europian allele counts
                case {'African', '_African'}
                    tmp_allele_counts = str2nums(tmp_str{3}); % Take African allele counts
            end
            if(i == 1) % get list of SNPs from first SNP
                field_ind_vec = zeros(length(S.field_names), 1); equal_sign_ind_vec = length_cell(S.field_names)+2;
                for j=1:length(S.field_names)
                    field_ind_vec(j) = strmatch([S.field_names{j} '='], tmp_str); % here we assume this never changes !!!!
                end
                gene_ind = field_ind_vec( strmatch('GL', S.field_names, 'exact') );
                feature_ind = field_ind_vec( strmatch('FG', S.field_names, 'exact') );
                
            end % if i == 1
            S.XXX_VARIANT_COUNT_(i) = tmp_allele_counts(1);
            S.XXX_REF_ALLELE_COUNT_(i) = tmp_allele_counts(2);
            S.XXX_FEATURE_{i} = str2word(',', tmp_str{feature_ind}(4:end), 1); % why this is hard-coded? 
            if(~isempty( strfind(S.XXX_FEATURE_{i}, ':') )) % remove gene name 
                S.XXX_FEATURE_{i} = str2word(':', S.XXX_FEATURE_{i}, 2); 
            end
            S.GENE{i} = tmp_str{gene_ind}(4:end); % get gene name
            for j=1:length(S.field_names) % loop on all fields - this is quite slow !! ('eval' inside a loop over all SNPs and all fields
                %%%%                field_ind = strmatch([S.field_names{j} '='], tmp_str);
                %%%%                eval_str = ['S.' S.field_descriptions{j} '{' num2str(i) '} = ''' ...
                %%%%                    str2word('=', tmp_str{field_ind_vec(j)}, 2) ''';'];
                S.INFO_ARR{i,j} = tmp_str{field_ind_vec(j)}(equal_sign_ind_vec(j):end);
                %%%%                eval_str = ['S.' S.field_descriptions{j} '{' num2str(i) '} = ''' ...
                %%%%                    tmp_str{field_ind_vec(j)}(equal_sign_ind_vec(j):end) ''';'];
                %%%%                eval(eval_str);
            end
        end % loop on SNPs
        S = my_rmfield(S, {'European', 'African', 'dbSNP', 'Total', 'Minor', ...
            'Observed', 'Average', 'geneList', 'Whether', 'PubMed', 'GenotypeMat', 'INFO'});  % why are these fields removed? save space!!!! (don't keep info field!!!)
        rm_field_names = {'EXOME_CHIP', 'GWAS_PUBMED'}; % in the future use this to remove unnecessary fields  
        if(isfield(S, 'scorePhastCons') && iscell(S.scorePhastCons))  % convert to .mat to save space and ease work
            S.scorePhastCons = cell2mat(empty_cell_to_numeric_val(str2num_cell(S.scorePhastCons), -9999999));
        end
        if(isfield(S, 'consScoreGERP') && iscell(S.consScoreGERP))
            S.consScoreGERP = cell2mat(empty_cell_to_numeric_val(str2num_cell(S.consScoreGERP), -9999999));
        end
        if(isfield(S, 'clinicalAssociation'))
            S.clinicalAssociation = strrep_cell(S.clinicalAssociation, 'unknown', '');
        end
        save([remove_suffix_from_file_name(site_frequency_file_name) population{1} '.mat'], '-struct', 'S'); % add new fields% remove fields to reduce memory
        return_flag=1 % only convert to .mat
    end % loop on population
    % % % % else % here .mat file already exists
    % % % %     for population = populations_vec
    % % % %         S = load([remove_suffix_from_file_name(site_frequency_file_name) population{1} '.mat'], ...
    % % % %             'XXX_VARIANT_COUNT_', 'XXX_REF_ALLELE_COUNT_', 'XXX_FEATURE_', 'XXX_GENE_', 'GENE');
    % % % %
    % % % %         if(exist('gene_list', 'var') && ~isempty(gene_list)) % temp hack to extract only a few genes
    % % % %             [~, keep_inds, J] = intersect_all(S.GENE, gene_list);
    % % % %             S = struct_by_inds(S, keep_inds);
    % % % %             save([remove_suffix_from_file_name(site_frequency_file_name) population{1} '.small_gene_list.mat'], '-struct', 'S');
    % % % %         end
    % % % %     end % loop on population
    %end % if .mat file exists
    if(return_flag) % leave function
        return;
    end
end % if extract_fields_flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% End Stage 2: Extract All Info Fields to Population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Stage 3: Compute Gene-Specific Matrices            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(compute_gene_matrices_flag)
    for population = populations_vec % perform further preprocessing (compute SNP specific parameters)
        S = load([remove_suffix_from_file_name(site_frequency_file_name) population{1} '.mat'], ... % load only neccessary fields
            'XXX_VARIANT_COUNT_', 'XXX_REF_ALLELE_COUNT_', 'XXX_FEATURE_', 'XXX_GENE_', 'GENE', ...
            'XXX_CHROM', 'POS'); % enable unique identifier for each allele
        
        % allele_types = {'Synonymous', 'NonSynonymous', 'intron', 'utr'};
        S.ALLELE_FREQ = S.XXX_VARIANT_COUNT_ ./ (S.XXX_VARIANT_COUNT_ + S.XXX_REF_ALLELE_COUNT_);
        [S.unique_genes, unique_gene_inds, S.GENE_INDS] = unique(S.GENE); S.num_genes = length(S.unique_genes);
        S.unique_chr = S.XXX_CHROM(unique_gene_inds);
        S.allele_types = unique(S.XXX_FEATURE_);
        S.num_allele_types = length(S.allele_types);
        n_vec = cell(1, S.num_allele_types); count_vec = cell(1, S.num_allele_types);
        f_vec = cell(1, S.num_allele_types); gene_inds = cell(1, S.num_allele_types);
        upper_freq_vec = [0.001 0.005 0.01 0.05 0.1 1]; % maximum allele freqeuncy
        
        if(compute_frac_carriers)
            T = cell(length(upper_freq_vec), 1);  allele_freq_vec = cell(length(upper_freq_vec), 1);
            for j=1:length(upper_freq_vec) % New: Compute frac. of carriers for each gene
                allele_freq_vec{j} = S.ALLELE_FREQ;
                allele_freq_vec{j}(allele_freq_vec{j} > upper_freq_vec(j)) = 0; % take only allele freq. below 1%
                T{j}.gene = S.unique_genes; % Save in tab-delimited format
            end
        end % compute frac. carriers
        S.num_alleles_per_gene_mat = zeros(S.num_allele_types, S.num_genes);
        S.total_heterozygosity_per_gene_mat = zeros(S.num_allele_types, S.num_genes);
        S.gene_by_allele_type_freq_list = cell(S.num_allele_types, S.num_genes); % NEW! list all allele frequencies in an array for each gene and each allele type
        S.gene_by_allele_type_n_list = cell(S.num_allele_types, S.num_genes); % NEW! list all number of alleles sequencied in an array for each gene and each allele type
        S.gene_by_allele_type_het_list = cell(S.num_allele_types, S.num_genes); % NEW! list all heterozygosities in an array for each gene and each allele type
        S.gene_by_allele_type_pos_list = cell(S.num_allele_types, S.num_genes); % NEW! Save identifiers (positions)
        
        S.gene_by_allele_type_inds_list = cell(S.num_allele_types, S.num_genes); % indices in original data
        
        for i=1:S.num_allele_types % Divide alleles to types
            sprintf('Create tables for allele type %ld out of %ld', i, S.num_allele_types)
            allele_type_inds = strfind_cell(lower(S.XXX_FEATURE_), lower(S.allele_types{i})); % find current alleles
            n_vec{i} =  S.XXX_REF_ALLELE_COUNT_(allele_type_inds) +  S.XXX_VARIANT_COUNT_(allele_type_inds);
            count_vec{i} =   S.XXX_VARIANT_COUNT_(allele_type_inds);
            f_vec{i} = S.XXX_VARIANT_COUNT_(allele_type_inds) ./ n_vec{i};
            
            [tmp_genes, ~, tmp_J] = unique(S.GENE(allele_type_inds) ); % get gene indices
            [~, ~, gene_inds{i}] = intersect(tmp_genes, S.unique_genes); % S.GENE(allele_type_inds), S.unique_genes);
            S.num_alleles_per_gene_mat(i,gene_inds{i}) = accumarray(tmp_J, ones(length(allele_type_inds), 1)); % get allele counts
            S.total_heterozygosity_per_gene_mat(i,gene_inds{i}) = accumarray(tmp_J, ...
                2.*S.ALLELE_FREQ(allele_type_inds) .* (1-S.ALLELE_FREQ(allele_type_inds))); % get total heterozygosity
            for j=1:length(gene_inds{i}) % loop on genes having this allele type
                S.gene_by_allele_type_inds_list{i, gene_inds{i}(j)} = allele_type_inds(tmp_J == j); % indices in original data
                S.gene_by_allele_type_freq_list{i, gene_inds{i}(j)} = S.ALLELE_FREQ(allele_type_inds(tmp_J == j));
                S.gene_by_allele_type_het_list{i, gene_inds{i}(j)} = 2 .* ...
                    S.gene_by_allele_type_freq_list{i, gene_inds{i}(j)} .* ...
                    (1 - S.gene_by_allele_type_freq_list{i, gene_inds{i}(j)});
                S.gene_by_allele_type_n_list{i, gene_inds{i}(j)} = S.XXX_VARIANT_COUNT_(allele_type_inds(tmp_J == j)) + ...
                    S.XXX_REF_ALLELE_COUNT_(allele_type_inds(tmp_J == j));
                S.gene_by_allele_type_pos_list{i, gene_inds{i}(j)} = S.POS(allele_type_inds(tmp_J == j)); % save positions
            end
            
            if(compute_frac_carriers)
                for j=1:length(upper_freq_vec) % loop on different frequency cutoffs
                    sum_freq_vec{i} = accumarray(tmp_J, allele_freq_vec{j}(allele_type_inds));
                    eval_str = ['T{' num2str(j)  '}.' strrep(S.allele_types{i}, '-', '_') ' = zeros(' num2str(S.num_genes) ', 1);']; eval(eval_str);
                    eval_str = ['T{' num2str(j)  '}.' strrep(S.allele_types{i}, '-', '_') '(gene_inds{' num2str(i) '}) = sum_freq_vec{' num2str(i) '};']; eval(eval_str);
                end
            end % if compute frac. carriers
        end % loop on allele types
        
        if(compute_frac_carriers) % get AGGREGATE allele types (containing multiple 'atomic' allele types) % if compute frac. carriers
            S.num_alleles_per_gene_mat = []; S.total_freq_per_gene_mat = []; S.total_heterozygosity_per_gene_mat = [];
            S.upper_freq_vec = upper_freq_vec;
            for j=1:length(upper_freq_vec) % compute freq. below threshold for different thresholds and classes
                sprintf('Compute frac. carriers for threshold %ld out of %ld', j, length(upper_freq_vec))
                [S.num_alleles_per_gene_mat{j} S.total_freq_per_gene_mat{j} ...
                    S.total_heterozygosity_per_gene_mat{j} ] = ...
                    get_cumulative_freq_internal(S, S.upper_freq_vec(j)); % this gives correct cumulative but doesn't collapse alleles !
                num_coding_allele_types = length(coding_allele_types);
                S.all_allele_types = S.allele_types;
                for i=1:num_coding_allele_types
                    eval_str = ['T{' num2str(j)  '}.' coding_allele_types{i}{1}  ' = zeros(' num2str(S.num_genes) ', 1);']; eval(eval_str);
                    
                    for k=1:length(coding_allele_types{i}{2})
                        if(~isfield (T{j}, coding_allele_types{i}{2}{k})) % some fields might be missing in certain chromosomes
                            eval_str = ['T{' num2str(j)  '}.' coding_allele_types{i}{2}{k}  ' = zeros(' num2str(S.num_genes) ', 1);']; eval(eval_str);
                        end
                        eval_str = ['T{' num2str(j)  '}.' coding_allele_types{i}{1} ' = T{' num2str(j)  '}.' coding_allele_types{i}{1} ...
                            ' + T{' num2str(j)  '}.' coding_allele_types{i}{2}{k} ';']; eval(eval_str);
                    end
                    [~, ~, J] = intersect(strrep_cell(S.allele_types, '-', '_'), coding_allele_types{i}{2});                   % Apply the same to S
                    S.num_alleles_per_gene_mat{j}(i+S.num_allele_types,:) = sum(S.num_alleles_per_gene_mat{j}(J,:));
                    S.total_freq_per_gene_mat{j}(i+S.num_allele_types,:) = sum(S.total_freq_per_gene_mat{j}(J,:));
                    S.total_heterozygosity_per_gene_mat{j}(i+S.num_allele_types,:) = sum(S.total_heterozygosity_per_gene_mat{j}(J,:));
                    
                    S.all_allele_types = [vec2row(S.all_allele_types) coding_allele_types{i}{1}];
                end
                S.num_all_allele_types = length(S.all_allele_types); % get #alleles with aggregate alleles
                output_filename = [remove_suffix_from_file_name(site_frequency_file_name) ...
                    population{1} '.freq_rare_variants_per_gene_DAF_less_' ...
                    strrep(num2str(upper_freq_vec(j)), '.', '_') ];
                WriteDataFile(T{j}  , [output_filename '.txt']);  % save in .txt format
                if(j == length(upper_freq_vec))
                    save([output_filename '.mat'], 'T'); % save in .mat format
                end
            end % loop on frequency threshold
        end % if compute frac carriers
        save([remove_suffix_from_file_name(site_frequency_file_name) population{1} '.mat'], '-append', '-struct', 'S'); % add new fields
        save([remove_suffix_from_file_name(site_frequency_file_name) population{1} '.mat'], '-append', ...
            'n_vec', 'count_vec', 'f_vec'); % , 'allele_types'); % Save again, add new fields
    end % loop on population again
end % if compute_gene_matrices_flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% End Stage 3: Compute Gene-Specific Matrices          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isempty(S)) % load data
    S = load( [remove_suffix_from_file_name(site_frequency_file_name) populations_vec{1} '.mat']); 
end


% Internal function: Get cumulative freq. of all alleles below a certain threshold
%
% Input:
% S - structure with all alleles
% freq_threshold - threshold below which to take alleles
%
% Output:
% num_alleles_per_gene_mat - table with total number of distinct alleles below threshold for each gene
% total_freq_per_gene_mat - table with cumulative frequency alleles below threshold for each gene
% total_heterozygosity_per_gene_mat - table with total heterozygosity of alleles below threshold for each gene
%
function [num_alleles_per_gene_mat, total_freq_per_gene_mat, total_heterozygosity_per_gene_mat] = ...
    get_cumulative_freq_internal(S, freq_threshold)

num_alleles_per_gene_mat = zeros(S.num_allele_types, S.num_genes);
total_freq_per_gene_mat =  zeros(S.num_allele_types, S.num_genes);
total_heterozygosity_per_gene_mat =  zeros(S.num_allele_types, S.num_genes);

for i=1:S.num_allele_types
    for j=1:S.num_genes
        I = find(S.gene_by_allele_type_freq_list{i,j} <= freq_threshold);
        num_alleles_per_gene_mat(i,j) = length(S.gene_by_allele_type_freq_list{i,j}(I));
        total_freq_per_gene_mat(i,j) = sum(S.gene_by_allele_type_freq_list{i,j}(I));
        total_heterozygosity_per_gene_mat(i,j) = 2.*sum( S.gene_by_allele_type_freq_list{i,j}(I) .* ...
            (1-S.gene_by_allele_type_freq_list{i,j}(I)) );
    end
end



