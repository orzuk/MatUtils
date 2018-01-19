% Master script for parsing data for site-frequency spectrum from multiple
% datasets and plot allele-frequencies and other results
if(isdeployed) % run as executable
    cd('/cs/cbio/orzuk/software/Code/Matlab');
    SetPathScript; % for running as executable
end
Assign24MammalsGlobalConstants; AssignGeneralConstants; AssignStatsConstants; AssignRVASConstants;
num_bins = 0:0.01:1; % bins for what?
params = ReadParametersFile(fullfile(github_dir, 'MatUtils', 'genetic_architecture', 'data', 'ExomeParsingParameters.txt'));
s_vec = [0 -logspace(-6, -1, 11)]; % set s-values for plotting and fitting
queue_str = 'priority'; % for submitting jobs at broad farm
global cumsum_log_vec;
cumsum_log_vec = cumsum([0 log(1:2*10000)]); % compute log-binomial coefficients to save time
%'rate.matrix.intergenic'; % matrix with codon annotations and mutation rates from Paz
exome_struct = get_exome_data_info(params.exome_data); % get metadata: file names, directories etc.

% for i=4:4 % Loop on datasets. Take only ESP data % length(spectrum_data_files)

%%%%%%%%%%%%%%%%%%%
% Parse SFS Files %
%%%%%%%%%%%%%%%%%%%
vcf_file_names =  GetFileNames(fullfile(spectrum_data_dir, exome_struct.sub_dir_str, [exome_struct.prefix, '*.vcf']), 1);
if(params.parse_site_frequency_flag) % here we parse
    %    if(read_to_mat_flag)
    for i=1:length(vcf_file_names) % 10 % TEMP!!! RUN ON FIRST 10 FILES FOR DEBUG. length(vcf_file_names) % loop on all chunks (By chromosomes or otherwise)
        parse_job_str = ['[A] =' ... % , n_vec, count_vec, f_vec, allele_types] = ' ...
            'parse_site_frequency_data(''' vcf_file_names{i} ...
            ''', exome_struct, [], ' num2str(params.read_to_mat_flag) ', ' num2str(params.extract_fields_flag) ', ' ...
            num2str(params.compute_gene_matrices_flag) ');']; %, gene_list
        if(in_matlab_flag)
            eval(parse_job_str);
        else
            SubmitMatlabJobToFarm(parse_job_str, ...
                fullfile('out', ['parse_' exome_struct.prefix '.' i '.out']), queue_str, [], [], mem_flag); % allow specifying memory allocation
        end
    end
end

%%%%%%%%%%%%%%%%%%
% Plot SFS Files %
%%%%%%%%%%%%%%%%%%
demography_file = [remove_suffix_from_file_name(exons_file) ...
    '_' 'AllPop' '_Demography.mat'];
demography_file = fullfile(spectrum_data_dir, ...
    exome_struct.data_str, 'AllPop', demography_file);
if(params.plot_site_frequency_flag) % Here we plot SFS for DATA !! for all populations !!!
    plot_site_frequency_data(fullfile(spectrum_data_dir, exome_struct.data_str, ...
        [exome_struct.prefix, '*.mat']), ... %  exome_struct.spectrum_data_file)  '.mat']), ... % '_' population{1}% _unique
        fullfile(mammals_data_dir, genome_version, [remove_suffix_from_file_name(exons_file) '_unique.mat']),  ... % GeneStruct
        exome_struct, ... %   {'European', 'African'}, ...
        fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file_exons), ... % TEMP! take old mutation rates file ! 
        [], [], [], [], exome_struct.target_length, num_bins, ...
        fullfile(spectrum_data_dir, 'out', exome_struct.data_str, exome_struct.prefix)); %   remove_suffix_from_file_name(exome_struct.spectrum_data_file)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Synonymous SNPs to fit demographic model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(params.fit_demography)
    if(~exist(fullfile(spectrum_data_dir, exome_struct.sub_dir_str, 'AllPop', [exome_struct.prefix '_AllPop.mat']), 'file'))
        for i=1:length(vcf_file_names) % loop on files to get a distribution from all chunks !!!
            iii = i
            spectrum_population_data_file = fullfile(dir_from_file_name(vcf_file_names{i}), ...
                'AllPop', [remove_suffix_from_file_name(remove_dir_from_file_name(vcf_file_names{i})) '_AllPop.mat']);
            A = load(spectrum_population_data_file, 'count_vec', 'f_vec', 'n_vec', 'allele_types', 'num_allele_types');
            if(i == 1)
                all_A = A;
            else
                all_A = union_SFS_structs(all_A, A);  % currently just unite all variants together - better to keep the gene identities information too? 
            end
        end
        save(fullfile(dir_from_file_name(spectrum_population_data_file), [exome_struct.prefix '_AllPop_union.mat']), '-struct', 'all_A');
    else
        all_A = load(fullfile(spectrum_data_dir, exome_struct.sub_dir_str, 'AllPop', [exome_struct.prefix '_AllPop_union.mat']));
    end
    all_A.mu = mu_per_site * 3*10^9 * 0.015 * 0.01 / 3; % TEMP!! estimated total mutation rate: mu_per_site * gene size / 3  for synonymous
    all_A.mu = all_A.mu * 1.5; % TEMP CORRECTION !!!
    synonymous_ind = find(strcmp( 'synonymous_variant', all_A.allele_types)); % 'synonymous_variant' % 'coding-synonymous'
    
    i_pop=0;
    for population = exome_struct.populations %  {'African'} % , 'African'} % European'} % ,
        i_pop=i_pop+1;
        fit_demography = 0;
        if(~exist(demography_file, 'file') || (1 == 0)) % here fit new model (computationally heavy)
            fit_demography = 1;
        else  % load model
            Demographic_model = cell(length(exome_struct.populations), 1); % temp allocate space
            load(demography_file);
            if((length(Demographic_model) < i_pop) || isempty(Demographic_model{i_pop}))
                fit_demography=1;
            end
            Demographic_model{i_pop}.name = ['Fitted.' population{1}];
        end
        if(fit_demography)
%            load(fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file)); % load mutation rates 
%            Unique = load(fullfile(mammals_data_dir, genome_version, [remove_suffix_from_file_name(exons_file) '_unique.mat'])); 
            % Alternative: read mutation rates from Samocha's paper: 
            MutationRateTable2 = load(fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file2));
            all_A.mu = sum(10.^(MutationRateTable2.syn)); %             all_A.mu = sum(UniqueMutationRateTable); % compute total mutation rate over ALL genes !! 
%            [GG, I_G, J_G] = intersect(upper(M.gene), upper(Unique.GeneStruct.gene_names));
%            figure; plot(10.^(M.syn(I_G)), MMM.UniqueMutationRateTable(J_G, 1), '*');
            
            [Demographic_model{i_pop}, max_LL_demographic_model(i_pop)] = ...
                fit_demographic_parameters_from_allele_spectrum( ...
                all_A.count_vec{synonymous_ind}(:,i_pop), all_A.n_vec{synonymous_ind}(:,i_pop), [],  all_A.mu); % all_A.mu(MutationTypes == SYNONYMOUS)); % fit demography
            Demographic_model{i_pop}.name = ['Fitted.' population{1}];
            save(demography_file, 'Demographic_model', 'max_LL_demographic_model'); % Save and plot demography
        end
        
        % Here plot all demographies:
        % 1. Plot the population size as function of generations for each demography
        % 2. Plot the SFS for different values of selection coefficient s for each demography
        
        load(demography_file);
        if(~isfield(Demographic_model{i_pop}, 'SFS') || (1 == 0)) %  ~isreal(Demographic_model{i_pop}.SFS.p_vec)) % add SFS to demographic mode
            %    s_vec = [0 -logspace(-6, -2, 4)]; % light run - just for debugging
            Demographic_model{i_pop}.iters = 10000; % number of alleles to simulate !!
            Demographic_model{i_pop}.s_grid = [0 -logspace(-6, -2, 101)]; % s vector for interpolation
            compute_flag = []; compute_flag.method = 'simulation'; compute_flag.smooth = 1;  Demographic_model{i_pop}.cond_on_polymorphic_flag=1
            [Demographic_model{i_pop}.SFS.x_vec, Demographic_model{i_pop}.SFS.p_vec, ...
                Demographic_model{i_pop}.SFS.L, SFS_compute_time] = ...
                compute_allele_freq_spectrum_from_demographic_model( ...
                Demographic_model{i_pop}, s_vec, compute_flag); %  s_vec([1 5 10])
            save(demography_file, 'Demographic_model', 'max_LL_demographic_model');
        end
        index_vec(i_pop) = Demographic_model{i_pop}.index;
    end % loop on populations (temp.)
    save(demography_file, '-append', 'index_vec'); % Save and plot demography
end % if fit demographies

if(params.plot_demography)
    % 1. Plot pop. size as function of time V
    % 2. plot SFS
    if(~exist('Demographic_model', 'var'))
         load(demography_file);
         all_A = load(fullfile(spectrum_data_dir, exome_struct.sub_dir_str, 'AllPop', [exome_struct.prefix '_AllPop_union.mat']), ...
             'n_vec', 'count_vec', 'allele_types');
        synonymous_ind = find(strcmp( 'synonymous_variant', all_A.allele_types)); % 'synonymous_variant' % 'coding-synonymous'
    end
    
    demographic_model_plot(Demographic_model, index_vec, max_LL_demographic_model, ...
        all_A.count_vec{synonymous_ind}, all_A.n_vec{synonymous_ind}, 0);  % plot properties of fitted demography: pop. size and SFS
end % if plot

% Temp: plot all populations together (should be part of plotting function
plot_params.figure_type = 1; plot_params.figs_dir = exome_data_figs_dir; plot_params.hist = 1; plot_params.xlim = [10^(-4) 1];
plot_params.cum=1; plot_params.weighted = 1; plot_params.normalize=1; plot_params.font_size=8; % plot cumulative weighted allele frequency distribution
if(~exist('Demographic_model', 'var'))
         load(demography_file);
end
plot_allele_freq(s_vec, Demographic_model(1:6), plot_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Missense+Stop SNPs to fit selection and tolerance parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit gene-specific selection and tolerance parameters             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(params.estimate_gene_by_gene) % estimate potential target size for each gene in the genome
    if(~exist(fullfile(mammals_data_dir, genome_version, exons_file), 'file')) % get all sequences
        GeneStruct = ExtractExons(mammals_data_dir, 'hg18', [], exons_file, 0); % Get gene sequences. (Don't get pwms!!!)
        save(fullfile(mammals_data_dir, genome_version, exons_file), '-struct', 'GeneStruct');
    else
        GeneStruct = load(fullfile(mammals_data_dir, genome_version, exons_file), ...
            'chr_vec', 'pos_start_vec', 'pos_end_vec', 'seqs', 'strand', 'gene_names', 'sort_perm'); % don't load sequences?
    end
    
    if(~exist(fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file), 'file')) % get estimated mutation rate per gene
        TripletsMutationTable = load(fullfile(spectrum_data_dir, 'mutation_rates', triplet_mutations_file)); % read 64x64 table
        [MutationRateTable, MutationTypes] = ComputeGeneMutationRates(TripletsMutationTable, ...
            fullfile(mammals_data_dir, genome_version, exons_file)); % Compute table of mutation rates for all genes
        save(fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file), 'MutationRateTable', 'MutationTypes');
    else
        %                load(fullfile(spectrum_data_dir, mutation_rates_file), 'MutationRateTable', 'MutationTypes');
    end
    % Need to loop here also on chunks !!
    for gene_prefix = {''} % {'ABCG1'} % for chrom 21 {'ANGPTL'} % for chrom 1 %%%% (num2cell(['A':'Z' '0':'9']'))'  %% {'ANKRD20A3'} %%  %% {'ANGP'} %% (mat2cell(['A':'Z' '0':'9']', ones(36,1), 1))' % enable also weird genes starting with a number
        gene_by_gene_job_str = ['parse_site_frequency_gene_by_gene(''' spectrum_data_dir ''', ''' exome_struct.spectrum_data_files_str ''', ' ... % spectrum_data_files{i}
            '''' fullfile(spectrum_data_dir, exome_data, 'GeneByGene') ''', ' ... % 'Tennessen_Science_2012'
            '''' fullfile(mammals_data_dir, genome_version, exons_file) ''' , ' ... % GeneStruct
            '''' fullfile(spectrum_data_dir, 'mutation_rates', mutation_rates_file) ''', ' ...
            '' '[], Demographic_model' ', ' ...
            num2str(params.plot_gene_by_gene) ', ' ...
            ' ''' gene_prefix{1} ''');'];
        if(in_matlab_flag)
            eval(gene_by_gene_job_str);
        else
            SubmitMatlabJobToFarm(gene_by_gene_job_str, ...
                fullfile('out', ['run_genes_prefix_' gene_prefix{1} '.out']), queue_str);
        end
    end % loop on prefix
end % estimate gene by gene parameters

% New! plot results of fitting : compare between populations, compare to human-chimp etc.
if(params.plot_gene_by_gene)
    plot_fitted_selection_parameters(fullfile(spectrum_data_dir, exome_data, 'GeneByGene', ...
        [remove_suffix_from_file_name(exons_file), '_fitted_stats.mat']), exome_struct, exome_data_figs_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % debug problems in generated SFS
% % % if(debug_sfs)
% % %     load('temp_surface.Fitted.European.mat');
% % %     SFS = load('temp_surface.Fitted.European.mat');
% % %     figure;     selection_color_vec = {'k', 'c', 'b', 'g', orange, 'r'}; % replace yellow with orange
% % %     my_symbol_vec = {'--', '-'}; % flip ordering (set integer powers as solid lines)
% % %     num_s = length(SFS.s_vec);
% % %     for i_s=1:length(SFS.s_vec)
% % %         tmp_color_ind = mod_max(ceil(i_s/2), ceil(num_s/2));
% % %         y_vec = cumsum_hist(SFS.x_vec, SFS.x_vec .* SFS.p_mat(i_s,:)); y_vec = y_vec ./ max(y_vec);
% % %         semilogx(SFS.x_vec ./ max(SFS.x_vec), y_vec, 'color', selection_color_vec{tmp_color_ind}, ...
% % %             'linestyle', my_symbol_vec{mod_max(i_s,2)}, 'linewidth', 2); hold on;
% % %         xlim([10^(-4) 1]);
% % %     end
% % %     legend(s_vec_to_legend(SFS.s_vec));
% % % end


% % Next two are already in previous gene-by-gene script?
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Aggregate s estimator for each gene from all populations %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(params.aggregate_population_estimators) % compute an aggregate esitmator for each gene from multiple populations
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Test for differences for s estimator for each gene from all populations %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(params.test_population_differences) % test for differences in selection between different populations for each gene
% end




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % Constants moved to parameters file:
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % exome_data = 'ExAC'; % 'ESP'; % 'ExAC'; % NEW! add also Exome Aggregation Data!!!!
% % % % % % % % % % % % % % % % Need to add also gnomad data !
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % old_run=0;
% % % % % % % % % % % % % % % % Set of flags determining which analysis to perform
% % % % % % % % % % % % % % % parse_site_frequency_flag = 0; % parse original datafile (different between different datasets)
% % % % % % % % % % % % % % % read_vcf_flag=0; % read vcf files for exome data
% % % % % % % % % % % % % % % unite_flag=0; % 0: parse ESP data. 1: unite all data to one chromosome
% % % % % % % % % % % % % % % read_to_mat_flag=0; % convert vcf (?) or other files to .mat format
% % % % % % % % % % % % % % % extract_fields_flag=1; % extract fields ???
% % % % % % % % % % % % % % % compute_gene_matrices_flag=1; % 1. Compute for each gene ?? flag for parsing ???
% % % % % % % % % % % % % % % plot_site_frequency_flag = 0; % 1: plot SFS data (this is also part of pre-processing)
% % % % % % % % % % % % % % % estimate_gene_by_gene = 1; % 1: analyze each gene seperately - estimate target size for each gene. This is what we want now!!!
% % % % % % % % % % % % % % % plot_gene_by_gene = 0; % make figures for individual genes
% % % % % % % % % % % % % % % fit_demography = 1;  % NEW! here fit a demographic model using only synonymous SNPs
% % % % % % % % % % % % % % % plot_demography = 1; % summary plots for demographic models
% % % % % % % % % % % % % % % aggregate_population_estimators = 0; % NEW! aggregate estimators from different populations
% % % % % % % % % % % % % % % test_population_differences = 0; % NEW! test for different in selection between different populations





