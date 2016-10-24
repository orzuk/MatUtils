% Master script for parsing data for site-frequency spectrum from multiple
% datasets and plot allele-frequencies and other results
Assign24MammalsGlobalConstants; AssignGeneralConstants; AssignStatsConstants; AssignRVASConstants;
num_bins = 0:0.01:1; % bins for what?

exome_data = 'ExAC'; % 'ESP'; % 'ExAC'; % NEW! add also Exome Aggregation Data!!!!

% Set of flags determining which analysis to perform
parse_site_frequency_flag = 1; % parse original datafile (different between different datasets)
read_vcf_flag=1; % read vcf files for exome data
unite_flag=0; % 0: parse ESP data. 1: unite all data to one chromosome
read_to_mat_flag=1; % convert vcf (?) or other files to .mat format
extract_fields_flag=0; % extract ??? fields
compute_gene_matrices_flag=0; % 1. Compute for each gene ?? flag for parsing ???
plot_site_frequency_flag = 1; % 1: plot SFS data (this is also part of pre-processing)
estimate_gene_by_gene = 0; % 1: analyze each gene seperately - estimate target size for each gene. This is what we want now!!!
plot_gene_by_gene = 0; % make figures for individual genes
fit_demography = 1;  % NEW! here fit a demographic model using only synonymous SNPs
aggregate_population_estimators = 0; % NEW! aggregate estimators from different populations
test_population_differences = 0; % NEW! test for different in selection between different populations

queue_str = 'priority'; % for submitting jobs at broad farm

global cumsum_log_vec;
cumsum_log_vec = cumsum([0 log(1:2*10000)]); % compute log-binomial coefficients to save time

%'rate.matrix.intergenic'; % matrix with codon annotations from Pazik

exome_struct = get_exome_data_info(exome_data); % get metadata: file names, directories etc.

% for i=4:4 % Loop on datasets. Take only ESP data % length(spectrum_data_files)
i_pop=1;


if(parse_site_frequency_flag) % here we parse 
    if(read_to_mat_flag)
        vcf_file_names =  GetFileNames(fullfile(spectrum_data_dir, exome_struct.sub_dir_str, [exome_struct.prefix, '*.vcf']), 1); 
        for i=1:length(vcf_file_names) % loop on all chunks (By chromosomes or otherwise)
              job_str = ['[A] =' ... % , n_vec, count_vec, f_vec, allele_types] = ' ...
                'parse_site_frequency_data(''' vcf_file_names{i} ...
                ''', [], ' num2str(read_to_mat_flag) ', ' num2str(extract_fields_flag) ', ' ...
                num2str(compute_gene_matrices_flag) ');']; %, gene_list

                if(in_matlab_flag)
                    eval(job_str);
                else
                    SubmitMatlabJobToFarm(job_str, ...
                        fullfile('out', ['parse_' exome_struct.prefix '.' i '.out']), queue_str, [], [], mem_flag); % allow specifying memory allocation
                end
        end
    end
end
        

for population = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
    if(~strcmp(population, 'African')) % temp: work only on one population1
        continue;
    end
    if(read_vcf_flag)
        max_chr=23; min_chr=1;
    else % run only once 
        max_chr=1; min_chr=1;
    end
    sub_dir_str = dir_from_file_name(exome_struct.spectrum_data_file);
    chr_file_str = suffix_from_file_name(remove_suffix_from_file_name(exome_struct.spectrum_data_file));
    
    %                min_chr=22; max_chr=23; % TEMP FOR DEBUG! TAKE SHORT CHROMOSOME
    for chr = min_chr:max_chr % Issue: not always good to split by chromosomes !! % take short chrom for debugging min_chr:max_chr % 1:23
        do_chr = chr
        if(read_vcf_flag) % One file per chromosome. (includes ALL populations)
            tmp_file_name = GetFileNames(fullfile(spectrum_data_dir, sub_dir_str, ...
                [exome_struct.prefix '*.chr' chr_num2str(chr) '.' chr_file_str '.vcf' ]));
            exome_struct.spectrum_data_file = fullfile(sub_dir_str, tmp_file_name{1});
            %                 spectrum_data_files{i} = fullfile(sub_dir_str, ...
            %                     ['ESP6500.chr' chr_num2str(chr) '.' chr_file_str '.vcf' ]); % '_' population{1} '.vcf']; %    .vcf'];
        else % One file. Already includes population string
            %                spectrum_data_files{i} = ['/Tennessen_Science_2012/all_chr_ESP6500.snps' '_' population{1} '.mat']; %    .vcf'];
            
            if(read_to_mat_flag)
                exome_struct.spectrum_data_file = fullfile(dir_from_file_name(exome_struct.spectrum_data_file), ...
                    ['all_chr_' exome_struct.prefix '.' chr_file_str '.vcf']); %  '_' population{1} '.mat']; %    .vcf'];
            else % already have a .mat file
                exome_struct.spectrum_data_file = fullfile(dir_from_file_name(exome_struct.spectrum_data_file), ...
                    ['all_chr_' exome_struct.prefix '.' chr_file_str '.mat']); %  '_' population{1} '.mat']; %    .vcf'];
            end
        end
        exome_struct.spectrum_data_files_str = [exome_struct.spectrum_data_files_str '''' ...
            remove_suffix_from_file_name(exome_struct.spectrum_data_file) '_' population{1} '.mat' ''','];
        
        if(parse_site_frequency_flag)
            job_str = ['[A] =' ... % , n_vec, count_vec, f_vec, allele_types] = ' ...
                'parse_site_frequency_data(''' fullfile(spectrum_data_dir, exome_struct.spectrum_data_file) ...
                ''', [], ' num2str(read_to_mat_flag) ', ' num2str(extract_fields_flag) ', ' ...
                num2str(compute_gene_matrices_flag) ');']; %, gene_list
            
            %                 [A n_vec count_vec f_vec allele_types] = ...
            %                     parse_site_frequency_data(fullfile(spectrum_data_dir, spectrum_data_files{i})); % , gene_list);
            
            if(unite_flag)
                %                    if(chr == min_chr)
                %                        A = load(fullfile(spectrum_data_dir, sub_dir_str, ...
                %                            ['all_chr_ESP6500.' chr_file_str '_'  population{1} '_up_to_chr' num2str(chr) '.mat'])); % load union
                %                    else % load current
                A = load([remove_suffix_from_file_name(fullfile(spectrum_data_dir, exome_struct.spectrum_data_file))  '_' population{1} '.mat']);
                A = my_rmfield(A, 'INFO_ARR');
                %                    end
                
                in_matlab_flag=1;
            else % don't unite chroms
                if(in_matlab_flag)
                    eval(job_str);
                else
                    SubmitMatlabJobToFarm(job_str, ...
                        fullfile('out', ['parse_ESP_chr' chr_num2str(chr) '.out']), queue_str, ...
                        [], [], mem_flag); % allow specifying memory allocation
                end
            end % if unite chroms
            if(read_vcf_flag && in_matlab_flag) % unite different files (even without 'unite chr')
                if(isfield(A, 'GENE'))
                    A.GENE = vec2column(A.GENE);
                end
                if(unite_flag)
                    if(chr == min_chr) % 1)
                        all_A = A;
                    else
                        field_names = fieldnames(A);
                        unite_field_names = intersect(fieldnames(A), {'XXX_VARIANT_COUNT_', 'XXX_REF_ALLELE_COUNT_', 'XXX_FEATURE_', 'GENE', 'XXX_CHROM', ...
                            'POS',  'ALLELE_FREQ',   'GENE_INDS', 'unique_genes'});
                        for j=1:length(unite_field_names)  % concatenate all chromosomes to one file (is this for population file? or one file for all populations?)
                            unite_str = ['all_A.' unite_field_names{j} ' = [all_A.' unite_field_names{j} ''' A.' unite_field_names{j} ''']'';'];
                            eval(unite_str)
                        end
                        [intercet_allele_types, I_types, J_types] = intersect(all_A.allele_types, A.allele_types); % Here unite cell-array. Problem! for different chromosomes might have different #allele_types and their encoding
                        for j=1:length(I_types)
                            all_A.n_vec{I_types(j)} = [all_A.n_vec{I_types(j)}' A.n_vec{J_types(j)}']';
                            all_A.f_vec{I_types(j)} = [all_A.f_vec{I_types(j)}' A.f_vec{J_types(j)}']';
                            all_A.count_vec{I_types(j)} = [all_A.count_vec{I_types(j)}' A.count_vec{J_types(j)}']';
                        end
                        [diff_allele_types, I_diff_types] = setdiff(A.allele_types, all_A.allele_types); % Here get different #allele_types and their encoding
                        for j=1:length(I_diff_types)% New alleles
                            all_A.n_vec{all_A.num_allele_types+j} = A.n_vec{I_diff_types(j)};
                            all_A.f_vec{all_A.num_allele_types+j} = A.f_vec{I_diff_types(j)};
                            all_A.count_vec{all_A.num_allele_types+j} = A.count_vec{I_diff_types(j)};
                            all_A.allele_types{all_A.num_allele_types+j} = A.allele_types{I_diff_types(j)}; % add allele types
                        end
                        all_A.num_allele_types = length(all_A.allele_types); % update # of allele types
                    end
                end % unite chr
                close all;
            end % if read vcf
            if(unite_flag)  % Here unite - why only first population? (Europeans?)
                save(fullfile(spectrum_data_dir, sub_dir_str, ...
                    ['all_chr_' exome_struct.prefix '.' chr_file_str '_'  population{1} '_up_to_chr' num2str(chr) '.mat']), '-struct', 'all_A'); % Save union
            end
        end % if parse SFS data
    end % loop on chr.
    
    if(plot_site_frequency_flag)
        plot_site_frequency_data(fullfile(spectrum_data_dir, ...
            [remove_suffix_from_file_name(exome_struct.spectrum_data_file)  '.mat']), ... % '_' population{1}% _unique
            fullfile(mammals_data_dir, genome_version, [remove_suffix_from_file_name(exons_file) '_unique.mat']),  ... % GeneStruct
            exome_struct.populations, ... %   {'European', 'African'}, ...
            fullfile(spectrum_data_dir, mutation_rates_file), ...
            [], [], [], [], exome_struct.target_length, num_bins, ...
            fullfile(spectrum_data_dir, 'out', remove_suffix_from_file_name(exome_struct.spectrum_data_file)));
        % %         plot_site_frequency_data(A, n_vec, count_vec, f_vec, allele_types, exome_struct.target_length, num_bins, ...
        % %             fullfile(spectrum_data_dir, 'out', remove_suffix_from_file_name(exome_struct.spectrum_data_file)));
    end
    i_pop=i_pop+1;
end % loop on population (temp.)


if(fit_demography) % here use Synonymous SNPs to fit demographic model
    spectrum_population_data_file{i_pop} = [remove_suffix_from_file_name(exome_struct.spectrum_data_file) '_' population{1} '.mat'];
    all_A = load(fullfile(spectrum_data_dir, spectrum_population_data_file{i_pop}), 'count_vec', 'f_vec', 'n_vec', 'allele_types');
    synonymous_ind = find(strcmp( 'coding-synonymous', all_A.allele_types))
    [Demographic_model, max_LL_demographic_model] = ...
        fit_demographic_parameters_from_allele_spectrum( ...
        all_A.count_vec{synonymous_ind}, all_A.n_vec{synonymous_ind});
end


spectrum_data_files_str = [spectrum_data_files_str(1:end-1) '}'];
if(estimate_gene_by_gene) % estimate potential target size for each gene in the genome
    if(~exist(fullfile(mammals_data_dir, genome_version, exons_file), 'file')) % get all sequences
        GeneStruct = ExtractExons(mammals_data_dir, 'hg18', [], exons_file, 0); % Get gene sequences. (Don't get pwms!!!)
        save(fullfile(mammals_data_dir, genome_version, exons_file), '-struct', 'GeneStruct');
    else
        GeneStruct = load(fullfile(mammals_data_dir, genome_version, exons_file), ...
            'chr_vec', 'pos_start_vec', 'pos_end_vec', 'seqs', 'strand', 'gene_names', 'sort_perm'); % don't load sequences?
    end
    
    
    if(~exist(fullfile(spectrum_data_dir, mutation_rates_file), 'file')) % get estimated mutation rate per gene
        TripletsMutationTable = load(fullfile(spectrum_data_dir, 'mutation_rates', triplet_mutations_file)); % read 64x64 table
        [MutationRateTable MutationTypes] = ComputeGeneMutationRates(TripletsMutationTable, ...
            fullfile(mammals_data_dir, genome_version, exons_file)); % Compute table of mutation rates for all genes
        save(fullfile(spectrum_data_dir, mutation_rates_file), 'MutationRateTable', 'MutationTypes');
    else
        %                load(fullfile(spectrum_data_dir, mutation_rates_file), 'MutationRateTable', 'MutationTypes');
    end
    
    %            for gene_prefix = {'APOA5'} % mat2cell(['A':'Z' '0':'9']', ones(36,1), 1) % enable also weird genes starting with a number
    
    for gene_prefix = {''} % {'ABCG1'} % for chrom 21 {'ANGPTL'} % for chrom 1 %%%% (num2cell(['A':'Z' '0':'9']'))'  %% {'ANKRD20A3'} %%  %% {'ANGP'} %% (mat2cell(['A':'Z' '0':'9']', ones(36,1), 1))' % enable also weird genes starting with a number
        job_str = ['parse_site_frequency_gene_by_gene(''' spectrum_data_dir ''', ' spectrum_data_files_str ', ' ... % spectrum_data_files{i}
            '''' fullfile(spectrum_data_dir, exome_data, 'GeneByGene') ''', ' ... % 'Tennessen_Science_2012'
            '''' fullfile(mammals_data_dir, genome_version, exons_file) ''' , ' ... % GeneStruct
            '''' fullfile(spectrum_data_dir, mutation_rates_file) ''', ' ...
            num2str(plot_gene_by_gene) ', ' ...
            '[], ''' gene_prefix{1} ''');'];
        
        %                 parse_site_frequency_gene_by_gene(spectrum_data_dir, spectrum_data_files{i}, ...
        %                     fullfile(spectrum_data_dir, 'Tennessen_Science_2012', 'GeneByGene'), ...  % assume ESP data
        %                     fullfile(mammals_data_dir, genome_version, exons_file), ... % GeneStruct
        %                     fullfile(spectrum_data_dir, mutation_rates_file), [], gene_prefix); % Estimate s and alpha for each gene in the genome - how???
        if(in_matlab_flag)
            eval(job_str);
        else
            SubmitMatlabJobToFarm(job_str, ...
                fullfile('out', ['run_genes_prefix_' gene_prefix{1} '.out']), queue_str);
            %                        fullfile(spectrum_data_dir, 'Tennessen_Science_2012', 'GeneByGene', 'out', ['run_genes_prefix_' gene_prefix '.out']), ...
            %                        queue_str);
        end
        
    end % loop on prefix
    
    if(aggregate_population_estimators) % compute an aggregate esitmator for each gene from multiple populations
        
        
    end
    
    if(test_population_differences) % test for differences in selection between different populations for each gene
        
        
    end
    
    
    
end % estimate gene by gene parameters


% end % loop on SFS data type (currently use ESP)









