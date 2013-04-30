% Look at all genomic regions, and determine conserved k-mers for each of them
% 
% The inputs to this function are:
%
% data_dir - where most of the data is (e.g. regions file)
% pwms_dir - where the blocks of pi matrices are held (from Manuel)
% regions_file - a file representing each region in the genome and its type
% omega_file - a file saying the omega conservation for each base in the genome
% pi_file - the file saying what is the pi value for each base in the genome
% chroms - what chromosomes should we run on
% k_vec - vector of size k values of kmers 
% conserved_output_file - here we save the conserved structure 
%
% Output: 
% ConservedKmersStruct - a structure containing all the conservation
% information for each kmer (counts, moments, pvalues)
%
% Note: Currently we do everything very crudely. There are many boundaries
% in the regions we use and we treat them either as consecutive regions or
% ignore them ...
function ConservedKmersStruct = GenomicRegionsKmersConservation(data_dir, pwms_dir, omega_dir, regions_file, ...
    chroms, k_vec, seq_flag, only_singletons_flag, conserved_output_file)

Assign24MammalsGlobalConstants; % call script for initilizing constants
cons_flag = USE_PI; % USE_OMEGA is bad and currently should not be used !!! 

genomic_regions = load(fullfile(data_dir, regions_file)); % load regions file

k_vec = union(k_vec, 1); % ALWAYS calculate singleton statistics (use for background for k-mers statistics)
num_bins = 1000; % currently irrelevant ...
ConservedKmersStruct = {};
for t_i = 1:length(genome_types)
    for k_i=1:length(k_vec) % finally loop over possible values of k
        k=k_vec(k_i);
        ConservedKmersStruct{t_i}.kmers_inds{k} = [];
        ConservedKmersStruct{t_i}.kmers_counts{k} = [];
        ConservedKmersStruct{t_i}.kmers_weights{k} = [];
        ConservedKmersStruct{t_i}.kmers_weights_sqr{k} = [];
        ConservedKmersStruct{t_i}.singletons_weights_cov{k} = zeros(4, 4, k);
        ConservedKmersStruct{t_i}.singletons_weights_counts{k} = zeros(4, 4, k);
    end
end

% Initilize the conservation structure
for t_i = 1:length(genome_types)
    ConservedKmersStruct{t_i}.NormalizationFlag = 0;
    for k_i=1:length(k_vec)
        k=k_vec(k_i)
        ConservedKmersStruct{t_i}.weights_total_count{k} = 0;
        ConservedKmersStruct{t_i}.weights_mean{k} = 0;
        ConservedKmersStruct{t_i}.weights_mean_sqr{k} = 0;
    end
end

for i=1:length(chroms) % loop over chromosomes
    chr = chroms(i)
    pi_files_list = GetFileNames(fullfile(pwms_dir, ['chr', chr_num2str(chroms(i))], 'mat', 'chr*.pi.mat')); % block pi .mat files
    %    omega_files_list = dir( fullfile(omega_dir, ['chr', chr_num2str(chroms(i))], 'mat', 'chr*.omegas.mat')); % block omega .mat files

    if(seq_flag == HUMAN) % load the human sequence
        genome_used = load(fullfile(data_dir, 'hg17', ['chr' chr_num2str(chr) '.mat']));
    end
        
    num_files = length(pi_files_list)
    for j=1:num_files  % divide the work into many blocks and then integrate them
        ttt = cputime;
        doing_file = j
        load(fullfile(pwms_dir, ['chr', chr_num2str(chroms(i))], 'mat', pi_files_list{j})); % load pi values
        if(cons_flag == USE_OMEGA)   % load matching omega values
            omega_file_name = [pi_files_list{j}(1:end-6) 'omegas.mat'];
            load(fullfile(omega_dir, ['chr', chr_num2str(chroms(i))], 'mat', omega_file_name)); 
        end
        if(~isempty(PWMS)) % Some files are simply empty
            if(cons_flag == USE_OMEGA) % CURRENTLY THIS FLAG IS NOT SUPPORTED !!! 
                [POS_INTER I J] = intersect(PWMS_POSITIONS, OMEGA_POSITIONS);         % perform intersection of positions
                weights = OMEGA(J); % take omega in the relevant positions
                seq = pwm_to_consensus(PWMS(I,:)'); % calculate the consensus
            else % here use log-odds of pi
                POS_INTER = PWMS_POSITIONS;
                weights = PWMS_LOGLIKE; 
                if(seq_flag == ANCESTRAL)
                    seq = pwm_to_consensus(PWMS'); % Currently: sequence is taken from PI consensus !!! could this bias results?
                else % take the human sequence ! could this bias results? 
                    seq = genome_used.Sequence(POS_INTER);
                end
            end
            tic; POS_GENOME_TYPES = StratifyToGenomicRegions(chr, POS_INTER, [], genome_version, genomic_regions); % toc; % check which genomic regions are contained in each one

            block.kmers_unique = {}; block.kmers_inds = {}; block.kmers_counts = {};
            block.kmers_weights = {}; block.kmers_weights_sqr = {}; block.kmers_weights_pvals = {};
            for t_i = 1:length(genome_types)
                region_type_inds = find(POS_GENOME_TYPES == t_i);
                if(~isempty(region_type_inds)) % for some blocks we may lack some regions ...
                    for k_i=1:length(k_vec) % finally loop over possible values of k
                        k=k_vec(k_i)
                        if(length(region_type_inds) >= k) % we cannot count smaller bins
                            tic; [block.kmers_unique{k} block.kmers_inds{k} block.kmers_counts{k} ...
                                block.kmers_weights{k} block.kmers_weights_sqr{k} block.kmers_weights_pvals{k} ...
                                block.weights_mean{k} block.weights_std{k} block.singletons_weights_cov{k}, ...
                                block.singletons_weights_counts{k}] = ...
                                KmersConservation(seq(region_type_inds), weights(region_type_inds), k, 1, only_singletons_flag);  % toc; % compute kmers statistics for this block
                            tic; ConservedKmersStruct = ...
                                UpdateBlockConservationStatistics(ConservedKmersStruct, block, t_i, k, only_singletons_flag); % toc % An auxillary function
                        end
                    end
                end
            end
        end
        ttt_one_file = cputime - ttt
    end % finished loop on files j
    
end

% Get the p-value for the whole blocks: this is not alway needed but whatever ... 
ConservedKmersStruct = ConservationCountsToPvalues(ConservedKmersStruct, genome_types, k_vec, only_singletons_flag, 'pvals'); 
save(fullfile(data_dir, 'hg17', conserved_output_file), 'ConservedKmersStruct');    % save conservation output struct

TTT = 9999; 



