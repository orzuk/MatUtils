% Parse data for site-frequency spectrum from 3 papers and plot allele-frequencies
Assign24MammalsGlobalConstants;
AssignGeneralConstants;
AssignStatsConstants; 
num_bins = 0:0.01:1;

parse_site_frequency_flag = 0; read_vcf=0; unite_chr=0; % parse ESP data
read_to_mat_flag=0; extract_fields_flag=0; compute_gene_matrices_flag=0; % flags for parsing

plot_site_frequency_flag = 0; % plot ESP data (this is also part of pre-processing)
estimate_gene_by_gene = 1; % analyze each gene seperately


queue_str = 'priority'; % for submitting jobs at broad farm

switch machine % Get directory
    case UNIX
        spectrum_data_dir = '/seq/orzuk/common_disease_model/data/SiteFrequencySpectra/';
        if(~exist('in_matlab_flag', 'var'))
            in_matlab_flag = 0;
        end
        if(~exist('mem_flag', 'var')) % || isempty(mem_flag))
            mem_flag = 8; % allow large memory
        end
    case PC
        spectrum_data_dir = 'C:\\research\common_disease_model\data\SiteFrequencySpectra';
        %        spectrum_data_dir = 'T:\\common_disease_model\data\SiteFrequencySpectra';
        in_matlab_flag = 1;
end
spectrum_data_files = {'/Nelson_Science_2012/Nelson_Science_data_table_S2.txt', ...
    '/Tennessen_Science_2012/all_chr_ESP6500.snps.vcf', ...
    '/Tennessen_Science_2012/all_chr_ESP6500.snps.small_gene_list.vcf', ...
    '/new_ESP/ESP6500SI*.snps_indels.vcf', ...
    []}; % Data file with exome sequencing results %  '/Tennessen_Science_2012/ESP6500.chr19.snps.vcf', []}; % S3
spectrum_labels = {'Nelson et al.', 'Tennessen et al.', 'Tennessen et al. few genes', 'Keinan et al.'};
target_length_vec = [ 50*10^6/100,  50*10^6, 50*10^6/100, 50*10^6]; % estimated # of nucleotides in target

exons_file = 'hg18_exons.mat'; % where to save exons
triplet_mutations_file = 'scone_hg17_m3_64x64_rates.txt'; % 'triplets_file_human_chimp_baboon.mat'; % file with mutation rates from all 64 codons (estimated from human-chimp-baboon alignment)
mutation_rates_file = 'human_genes_mutation_rates_hg18.mat'; % output file with triplet mutation rates

%'rate.matrix.intergenic'; % matrix with codon annotations from Pazik

spectrum_data_files_str = '{';
for i=4:4 % Loop on datasets. Take only ESP data % length(spectrum_data_files)
    for population = {'European'} % , 'African'} % European'} % , 
        if(read_vcf)
            max_chr=23; min_chr=1;
        else
            max_chr=1; min_chr=1;
        end
        sub_dir_str = dir_from_file_name(spectrum_data_files{i});
        chr_file_str = suffix_from_file_name(remove_suffix_from_file_name(spectrum_data_files{i}));
        
        for chr = min_chr:max_chr % 1:23
            do_chr = chr
            if(read_vcf) % One file per chromosome. (includes ALL populations)
                tmp_file_name = GetFileNames(fullfile(spectrum_data_dir, sub_dir_str, ['ESP6500*.chr' chr_num2str(chr) '.' chr_file_str '.vcf' ]));
                spectrum_data_files{i} = fullfile(sub_dir_str, tmp_file_name{1});
%                 spectrum_data_files{i} = fullfile(sub_dir_str, ...
%                     ['ESP6500.chr' chr_num2str(chr) '.' chr_file_str '.vcf' ]); % '_' population{1} '.vcf']; %    .vcf'];
            else % One file. Already includes population string
                %                spectrum_data_files{i} = ['/Tennessen_Science_2012/all_chr_ESP6500.snps' '_' population{1} '.mat']; %    .vcf'];
                spectrum_data_files{i} = fullfile(dir_from_file_name(spectrum_data_files{i}), ...
                    ['all_chr_ESP6500.' chr_file_str '.mat']); %  '_' population{1} '.mat']; %    .vcf'];
            end
            spectrum_data_files_str = [spectrum_data_files_str '''' ...
                remove_suffix_from_file_name(spectrum_data_files{i}) '_' population{1} '.mat' ''','];
            
            if(parse_site_frequency_flag)
                job_str = ['[A] =' ... % , n_vec, count_vec, f_vec, allele_types] = ' ...
                    'parse_site_frequency_data(''' fullfile(spectrum_data_dir, spectrum_data_files{i}) ...
                    ''', [], ' num2str(read_to_mat_flag) ', ' num2str(extract_fields_flag) ', ' ...
                    num2str(compute_gene_matrices_flag) ');']; %, gene_list
                
                %                 [A n_vec count_vec f_vec allele_types] = ...
                %                     parse_site_frequency_data(fullfile(spectrum_data_dir, spectrum_data_files{i})); % , gene_list);
                
                if(unite_chr)
%                     if(chr == min_chr)
%                         A = load(fullfile(spectrum_data_dir, sub_dir_str, ...
%                             ['all_chr_ESP6500.' chr_file_str '_'  population{1} '_up_to_chr' num2str(chr) '.mat'])); % load union
%                     else
                        A = load([remove_suffix_from_file_name(fullfile(spectrum_data_dir, spectrum_data_files{i}))  '_' population{1} '.mat']);
%                     end
                    
                    in_matlab_flag=1;
                else
                    if(in_matlab_flag)
                        eval(job_str);
                    else
                        SubmitMatlabJobToFarm(job_str, ...
                            fullfile('out', ['parse_ESP_chr' chr_num2str(chr) '.out']), queue_str, ...
                            [], [], mem_flag); % allow specifying memory allocation
                    end
                end
                if(read_vcf && in_matlab_flag) % unite different files
                    if(isfield(A, 'GENE'))
                        A.GENE = vec2column(A.GENE);
                    end
                    if(chr == min_chr) % 1)
                        all_A = A;
                    else
                        field_names = fieldnames(A);
                        for j=1:length(field_names)  % concatenate all chromosomes to one file
                            eval_str = ['all_A.' field_names{j} ' = [all_A.' field_names{j} ''' A.' field_names{j} ''']'';'];
                            eval(eval_str)
                        end
                    end
                    close all;
                end % if read vcf
                if(unite_chr)
                    save(fullfile(spectrum_data_dir, sub_dir_str, ...
                        ['all_chr_ESP6500.' chr_file_str '_'  population{1} '_up_to_chr' num2str(chr) '.mat']), '-struct', 'all_A'); % Save union
                end
            end % if parse SFS data
        end % loop on chr.
        
        if(plot_site_frequency_flag)
            plot_site_frequency_data(fullfile(spectrum_data_dir, ...
                [remove_suffix_from_file_name(spectrum_data_files{i})  '.mat']), ... % '_' population{1}% _unique
                fullfile(mammals_data_dir, genome_version, [remove_suffix_from_file_name(exons_file) '_unique.mat']),  ... % GeneStruct
                {'European', 'African'}, ...
                fullfile(spectrum_data_dir, mutation_rates_file), ...
                [], [], [], [], target_length_vec(i), num_bins, ...
                fullfile(spectrum_data_dir, 'out', remove_suffix_from_file_name(spectrum_data_files{i})));
            % %         plot_site_frequency_data(A, n_vec, count_vec, f_vec, allele_types, target_length_vec(i), num_bins, ...
            % %             fullfile(spectrum_data_dir, 'out', remove_suffix_from_file_name(spectrum_data_files{i})));
        end
        
    end % loop on population (temp.)
    
    spectrum_data_files_str = [spectrum_data_files_str(1:end-1) '}'];
    if(estimate_gene_by_gene) % estimate potential target size for each gene in the genome
        if(~exist(fullfile(mammals_data_dir, genome_version, exons_file), 'file')) % get all sequences
            GeneStruct = ExtractExons(mammals_data_dir, 'hg18', [], exons_file, 0); % Get gene sequences. (Don't get pwms!!!)
        else
            %                 GeneStruct = load(fullfile(mammals_data_dir, genome_version, exons_file), ...
            %                     'chr_vec', 'pos_start_vec', 'pos_end_vec', 'seqs', 'strand', 'gene_names', 'sort_perm');
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
        
        for gene_prefix = {'ANGPTL'} %%%% (num2cell(['A':'Z' '0':'9']'))'  %% {'ANKRD20A3'} %%  %% {'ANGP'} %% (mat2cell(['A':'Z' '0':'9']', ones(36,1), 1))' % enable also weird genes starting with a number
            job_str = ['parse_site_frequency_gene_by_gene(''' spectrum_data_dir ''', ' spectrum_data_files_str ', ' ... % spectrum_data_files{i}
                '''' fullfile(spectrum_data_dir, 'Tennessen_Science_2012', 'GeneByGene') ''', ' ...
                '''' fullfile(mammals_data_dir, genome_version, exons_file) ''' , ' ... % GeneStruct
                '''' fullfile(spectrum_data_dir, mutation_rates_file) ''', ' ...
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
    end % estimate gene by gene parameters
    
    
end % loop on SFS data type (currently use ESP)









