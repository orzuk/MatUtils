% Master script for parsing data for site-frequency spectrum from multiple
% datasets and plot allele-frequencies and other results
Assign24MammalsGlobalConstants; AssignGeneralConstants; AssignStatsConstants; AssignRVASConstants;
num_bins = 0:0.01:1; % bins for what?

exome_data = 'ExAC'; % 'ESP'; % 'ExAC'; % NEW! add also Exome Aggregation Data!!!!
% Need to add also gnomad data !

old_run=0;
% Set of flags determining which analysis to perform
parse_site_frequency_flag = 0; % parse original datafile (different between different datasets)
read_vcf_flag=0; % read vcf files for exome data
unite_flag=0; % 0: parse ESP data. 1: unite all data to one chromosome
read_to_mat_flag=0; % convert vcf (?) or other files to .mat format
extract_fields_flag=1; % extract ??? fields
compute_gene_matrices_flag=1; % 1. Compute for each gene ?? flag for parsing ???
plot_site_frequency_flag = 0; % 1: plot SFS data (this is also part of pre-processing)
estimate_gene_by_gene = 1; % 1: analyze each gene seperately - estimate target size for each gene. This is what we want now!!!
plot_gene_by_gene = 0; % make figures for individual genes
fit_demography = 0;  % NEW! here fit a demographic model using only synonymous SNPs
plot_demographies = 0; % summary plots for demographic models 
aggregate_population_estimators = 0; % NEW! aggregate estimators from different populations
test_population_differences = 0; % NEW! test for different in selection between different populations

queue_str = 'priority'; % for submitting jobs at broad farm
global cumsum_log_vec;
cumsum_log_vec = cumsum([0 log(1:2*10000)]); % compute log-binomial coefficients to save time
%'rate.matrix.intergenic'; % matrix with codon annotations and mutation rates from Paz
exome_struct = get_exome_data_info(exome_data); % get metadata: file names, directories etc.

% for i=4:4 % Loop on datasets. Take only ESP data % length(spectrum_data_files)
i_pop=1;


%%%%%%%%%%%%%%%%%%%
% Parse SFS Files %
%%%%%%%%%%%%%%%%%%%
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
end

% return;


%%%%%%%%%%%%%%%%%%
% Plot SFS Files %
%%%%%%%%%%%%%%%%%%
demography_file = [remove_suffix_from_file_name(exons_file) ...
    '_' 'AllPop' '_Demography.mat'];
demography_file = fullfile(spectrum_data_dir, ...
    exome_struct.data_str, 'AllPop', demography_file);
for population = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
    if(~strcmp(population, 'African')) % temp: work only on one population!
        continue;
    end
    if(plot_site_frequency_flag) % Here we plot SFS for DATA !! 
        plot_site_frequency_data(fullfile(spectrum_data_dir, exome_struct.data_str, ...
            [exome_struct.prefix, '*.mat']), ... %  exome_struct.spectrum_data_file)  '.mat']), ... % '_' population{1}% _unique
            fullfile(mammals_data_dir, genome_version, [remove_suffix_from_file_name(exons_file) '_unique.mat']),  ... % GeneStruct
            exome_struct, ... %   {'European', 'African'}, ...
            fullfile(spectrum_data_dir, mutation_rates_file), ...
            [], [], [], [], exome_struct.target_length, num_bins, ...
            fullfile(spectrum_data_dir, 'out', exome_struct.data_str, exome_struct.prefix)); %   remove_suffix_from_file_name(exome_struct.spectrum_data_file)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use Synonymous SNPs to fit demographic model %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(fit_demography)
        spectrum_population_data_file{i_pop} = fullfile(dir_from_file_name(exome_struct.spectrum_data_file), ...
            population{1}, [remove_suffix_from_file_name(remove_dir_from_file_name(exome_struct.spectrum_data_file)) '_' population{1} '.mat']);
        all_A = load(fullfile(spectrum_data_dir, spectrum_population_data_file{i_pop}), 'count_vec', 'f_vec', 'n_vec', 'allele_types');
        all_A.mu = mu_per_site * 3*10^9 * 0.015 * 0.01 / 3; % TEMP!! estimated total mutation rate: mu_per_site * gene size / 3  for synonymous
        all_A.mu = all_A.mu * 1.5; % TEMP CORRECTION !!!
        synonymous_ind = find(strcmp( 'synonymous_variant', all_A.allele_types)) % 'synonymous_variant' % 'coding-synonymous'
        if(~exist(demography_file, 'file'))
            [Demographic_model{i_pop}, max_LL_demographic_model(i_pop)] = ...
                fit_demographic_parameters_from_allele_spectrum( ...
                all_A.count_vec{synonymous_ind}, all_A.n_vec{synonymous_ind}, [],  all_A.mu); % PROBLEM HERE!! WORK (my implementation / software)
            Demographic_model{i_pop}.name = ['Fitted.' population{1}];
            save(demography_file, 'Demographic_model', 'max_LL_demographic_model'); % Save and plot demography
        else
            load(demography_file);
            Demographic_model{i_pop}.name = ['Fitted.' population{1}];
        end
        demographic_model_plot(Demographic_model, Demographic_model.index, max_LL_demographic_model, ...
            all_A.count_vec{synonymous_ind}, all_A.n_vec{synonymous_ind}, 0);
    end
    
    i_pop=i_pop+1;
end % loop on population (temp.)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Missense+Stop SNPs to fit selection and tolerance parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(demography_file);
spectrum_data_files_str = [exome_struct.spectrum_data_files_str(1:end-1) '}'];

if(~isfield(Demographic_model{i_pop}, 'SFS')) % add SFS to demographic mode
%    s_vec = [0 -logspace(-6, -2, 4)]; % light run - just for debugging 
    s_vec = [0 -logspace(-6, -1, 11)]; % heavier run 
    Demographic_model{i_pop}.iters = 50; 
    Demographic_model{i_pop}.s_grid = [0 -logspace(-6, -2, 101)]; % s vector for interpolation
    compute_flag = []; compute_flag.method = 'simulation'; compute_flag.smooth = 1;
    [Demographic_model{i_pop}.SFS.x_vec, Demographic_model{i_pop}.SFS.p_vec, ...
        Demographic_model{i_pop}.SFS.L, SFS_compute_time] = ...
        compute_allele_freq_spectrum_from_demographic_model( ...
        Demographic_model{i_pop}, s_vec, compute_flag); 
    
    save(demography_file, 'Demographic_model', 'max_LL_demographic_model');    
end

% Here plot all demographies: 
% 1. Plot the population size as function of generations for each demography
% 2. Plot the SFS for different values of selection coefficient s for each demography
if(plot_demographies) 
    % 1. Plot pop. size
    
    % 2. plot SFS
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit gene-specific selection and tolerance parameters             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        [MutationRateTable, MutationTypes] = ComputeGeneMutationRates(TripletsMutationTable, ...
            fullfile(mammals_data_dir, genome_version, exons_file)); % Compute table of mutation rates for all genes
        save(fullfile(spectrum_data_dir, mutation_rates_file), 'MutationRateTable', 'MutationTypes');
    else
        %                load(fullfile(spectrum_data_dir, mutation_rates_file), 'MutationRateTable', 'MutationTypes');
    end
    
    %            for gene_prefix = {'APOA5'} % mat2cell(['A':'Z' '0':'9']', ones(36,1), 1) % enable also weird genes starting with a number
    
    % Need to loop here also on chunks !! 
    for gene_prefix = {''} % {'ABCG1'} % for chrom 21 {'ANGPTL'} % for chrom 1 %%%% (num2cell(['A':'Z' '0':'9']'))'  %% {'ANKRD20A3'} %%  %% {'ANGP'} %% (mat2cell(['A':'Z' '0':'9']', ones(36,1), 1))' % enable also weird genes starting with a number
        job_str = ['parse_site_frequency_gene_by_gene(''' spectrum_data_dir ''', ' spectrum_data_files_str ', ' ... % spectrum_data_files{i}
            '''' fullfile(spectrum_data_dir, exome_data, 'GeneByGene') ''', ' ... % 'Tennessen_Science_2012'
            '''' fullfile(mammals_data_dir, genome_version, exons_file) ''' , ' ... % GeneStruct
            '''' fullfile(spectrum_data_dir, mutation_rates_file) ''', ' ...
            '' '[], Demographic_model' ', ' ...
            num2str(plot_gene_by_gene) ', ' ...
            ' ''' gene_prefix{1} ''');'];        
        %                 parse_site_frequency_gene_by_gene(spectrum_data_dir, spectrum_data_files{i}, ...
        %                     fullfile(spectrum_data_dir, 'Tennessen_Science_2012', 'GeneByGene'), ...  % assume ESP data
        %                     fullfile(mammals_data_dir, genome_version, exons_file), ... % GeneStruct
        %                     fullfile(spectrum_data_dir, mutation_rates_file), [], gene_prefix); % Estimate s and alpha for each gene in the genome - how???
        if(in_matlab_flag)
            eval(job_str);
        else
            SubmitMatlabJobToFarm(job_str, ...
                fullfile('out', ['run_genes_prefix_' gene_prefix{1} '.out']), queue_str);
        end
        
    end % loop on prefix
    
    
    % Next two are already in previous gene-by-gene script? 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aggregate s estimator for each gene from all populations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(aggregate_population_estimators) % compute an aggregate esitmator for each gene from multiple populations
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test for differences for s estimator for each gene from all populations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(test_population_differences) % test for differences in selection between different populations for each gene
     
    end
end % estimate gene by gene parameters


% end % loop on SFS data type (currently use ESP)









