% Master script for parsing data for site-frequency spectrum from multiple
% datasets and plot allele-frequencies and other results
Assign24MammalsGlobalConstants; AssignGeneralConstants; AssignStatsConstants; AssignRVASConstants;
num_bins = 0:0.01:1; % bins for what?

exome_data = 'ExAC'; % 'ESP'; % 'ExAC'; % NEW! add also Exome Aggregation Data!!!!

old_run=0;
% Set of flags determining which analysis to perform
parse_site_frequency_flag = 1; % parse original datafile (different between different datasets)
read_vcf_flag=0; % read vcf files for exome data
unite_flag=0; % 0: parse ESP data. 1: unite all data to one chromosome
read_to_mat_flag=1; % convert vcf (?) or other files to .mat format
extract_fields_flag=1; % extract ??? fields
compute_gene_matrices_flag=0; % 1. Compute for each gene ?? flag for parsing ???
plot_site_frequency_flag = 0; % 1: plot SFS data (this is also part of pre-processing)
estimate_gene_by_gene = 0; % 1: analyze each gene seperately - estimate target size for each gene. This is what we want now!!!
plot_gene_by_gene = 0; % make figures for individual genes
fit_demography = 0;  % NEW! here fit a demographic model using only synonymous SNPs
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
    %    if(read_to_mat_flag)
    vcf_file_names =  GetFileNames(fullfile(spectrum_data_dir, exome_struct.sub_dir_str, [exome_struct.prefix, '*.vcf']), 1);
    for i=1:10 % TEMP!!! RUN ON FIRST 10 FILES FOR DEBUG. length(vcf_file_names) % loop on all chunks (By chromosomes or otherwise)
        job_str = ['[A] =' ... % , n_vec, count_vec, f_vec, allele_types] = ' ...
            'parse_site_frequency_data(''' vcf_file_names{i} ...
            ''', exome_struct, [], ' num2str(read_to_mat_flag) ', ' num2str(extract_fields_flag) ', ' ...
            num2str(compute_gene_matrices_flag) ');']; %, gene_list
        
        if(in_matlab_flag)
            eval(job_str);
        else
            SubmitMatlabJobToFarm(job_str, ...
                fullfile('out', ['parse_' exome_struct.prefix '.' i '.out']), queue_str, [], [], mem_flag); % allow specifying memory allocation
        end
    end
    
    %    else % here
    
    %    end
end

% return;


for population = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
    
    if(~strcmp(population, 'African')) % temp: work only on one population!
        continue;
    end
    
    if(old_run)
        old_run_parse_site_frequency_data;
    end
    
    if(plot_site_frequency_flag)
        plot_site_frequency_data(fullfile(spectrum_data_dir, exome_struct.data_str, ...
            [exome_struct.prefix, '*.mat']), ... %  exome_struct.spectrum_data_file)  '.mat']), ... % '_' population{1}% _unique
            fullfile(mammals_data_dir, genome_version, [remove_suffix_from_file_name(exons_file) '_unique.mat']),  ... % GeneStruct
            exome_struct.populations, ... %   {'European', 'African'}, ...
            fullfile(spectrum_data_dir, mutation_rates_file), ...
            [], [], [], [], exome_struct.target_length, num_bins, ...
            fullfile(spectrum_data_dir, 'out', exome_struct.data_str, exome_struct.prefix)); %   remove_suffix_from_file_name(exome_struct.spectrum_data_file)));
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









