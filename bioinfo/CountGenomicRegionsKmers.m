% Count kmers in a set of coordinates given
% from the genome. Since the genome is big, the function does so
% by dividing into many blocks
% The prefix is used to count only kmers which share this prefix. This is
% good especially to save memory - we divide the genome into different
% files such that each file contains only a fraction of all kmers.
% We currently use prefix of length 4, so this divides the genome into
% 4^4=256 sets of kmers, each with roughly 10 million sequences
%
% We have two modes: One is regions which are assumed to be small,
% and the other is the whole genome (in case the pos_start_vec is empty)
%
% Input:
% chr_vec, pos_start_vec, pos_end_vec - specifying the regions (could be empty)
% k - The kmers size
% block_size - size of blocks we divide the genome to
% save_blocks_flag - 1 - count kmers and generate blocks. 0 - just collect blocks of kmers statistics and combine them
% prefix - collect only blocks with a given prefix (should be 0 if save_blocks_flag=1)
% conservation_flag - flag saying if also to perform conservation statistics
% data_dir - directory of input data
% organism_dir - subdirectory with the data, corresponding to genome version
% output_dir - where to write the output file/s
% chr_to_run - we may choose to run only specific chromosomes
% output_file - name of file to save results in
%
function [num_blocks regions_len kmers_unique kmers_counts kmers_weights kmers_weights_sqr] = ...
    CountGenomicRegionsKmers(chr_vec, pos_start_vec, pos_end_vec, k, block_size, save_blocks_flag, ...
    prefix, cumulative_flag, conservation_flag, ...
    data_dir, organism_dir, output_dir, chr_to_run, output_file)

genome_version = remove_dir_from_file_name(organism_dir) % extract the genome version here 

Assign24MammalsGlobalConstants;

if(isempty(chr_vec))
    chr_vec = [1:organism_chr+1];
end
chroms = unique(chr_vec); % see what chromosomes we have to consider
chroms = intersect(chroms, chr_to_run);
kmers_unique = []; kmers_counts = []; kmers_weights = []; kmers_weights_sqr = [];

chroms_is = chroms
organism_str_is = organism_str
output_dir_is = output_dir
output_file_is = output_file
conservation_is = conservation_flag
% kkk = rrr
%  return;

for i=1:length(chroms)
    if(cumulative_flag == 0) % in this case work on each chromosome seperately - faster and more parallelible but need to unite results at the end
        kmers_unique = []; kmers_counts = [];
    end

    chr = chroms(i);
    organism_is = organism_str;
    if(~exist(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)]), 'dir'))
        mkdir(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)]));
    end
    %%    if(save_blocks_flag == 1) % if not, we've already got the sequence splitted to blocks
    load(fullfile(data_dir, organism_dir, ['chr' chr_num2str(chr, organism_str) '.mat'])); % load chrom genome sequence file ..
    if(~exist('upper_flag', 'var'))  % just make sure everything is upper-case
        Sequence = upper(Sequence);
        upper_flag = 1; % if you want to keep it case sensitive, then upper_flag should be zero
        save(fullfile(data_dir, organism_dir, ['chr' chr_num2str(chr, organism_str) '.mat']), 'Sequence', 'upper_flag', '-append');
    end
    if(upper_flag == 0)
        Sequence = upper(Sequence);
    end
    %%    end

    if(~isempty(pos_start_vec)) % here we count the whole genome
        block_start_pos = pos_start_vec(chr_vec == chr);
        block_end_pos = pos_end_vec(chr_vec == chr);
    else  % here just use pre-defined regions
        num_blocks = ceil(length(Sequence) / block_size);
        block_start_pos = [0:num_blocks-1]*block_size+1;
        block_end_pos = [1:num_blocks]*block_size; block_end_pos(end) = min(block_end_pos(end), length(Sequence));
    end
    % throw away blocks which are shorter than k itself
    good_block_inds = find(block_end_pos - block_start_pos >= k); num_blocks = length(good_block_inds);
    block_start_pos = block_start_pos(good_block_inds); block_end_pos = block_end_pos(good_block_inds);
    regions_len = sum(block_end_pos - block_start_pos + 1);

    if((save_blocks_flag == 1) && conservation_flag) % see what files contain the blocks of genome pi values
        PI = ExtractPWMsByPositions(chr+zeros(1,length(good_block_inds)), block_start_pos, block_end_pos, [], genome_version, ...
            0, 1, 0, 0, 0); % save only the log-likelihood
        % see where we actually don't have conservation data ..
    end

    tic;  mid_block_kmers_unique = []; mid_block_kmers_counts = [];
    for j=1:num_blocks
        doing_block = j;
        if(save_blocks_flag == 1)% Just save the block and get out of here. Here output_file is not used
            if(conservation_flag) % here also compute the conservation
                positive_inds = find(PI.loglike{j} > 0);
                %                if(~isempty(positive_inds))
                if(length(positive_inds) >= k) % we require to have enough nucleotides which are alignable
                    block_seq = Sequence(PI.pos_start_vec(j):PI.pos_end_vec(j));
                    [block_kmers_unique block_kmers_inds block_kmers_counts block_kmers_weights block_kmers_weights_sqr] = ...
                        KmersConservation(block_seq(positive_inds), PI.loglike{j}(positive_inds), k, 0, 0);
                    save(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)], ...
                        ['block_kmers_' num2str(PI.pos_start_vec(j)) '_' num2str(PI.pos_end_vec(j)) '.mat']), ...
                        'block_kmers_unique', 'block_kmers_counts', 'block_kmers_weights', 'block_kmers_weights_sqr'); % just save the blocks
                end
            else
                [block_kmers_unique block_kmers_counts] = nmercount(Sequence(block_start_pos(j):block_end_pos(j)), k);
                save(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)], ...
                    ['block_kmers_' num2str(block_start_pos(j)) '_' num2str(block_end_pos(j)) '.mat']), ...
                    'block_kmers_unique', 'block_kmers_counts'); % just save the blocks
            end
        else  % read the blocks and unite them
            if(exist(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)], ...
                    ['block_kmers_' num2str(block_start_pos(j)) '_' num2str(block_end_pos(j)) '.mat']), 'file')) % just load the blocks
                load(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)], ...
                    ['block_kmers_' num2str(block_start_pos(j)) '_' num2str(block_end_pos(j)) '.mat'])); % just load the blocks
                num_block_kmers = size(block_kmers_unique, 1);
                [I J] = find( (block_kmers_unique == 'A') | (block_kmers_unique == 'C') | ...
                    (block_kmers_unique == 'G') | (block_kmers_unique == 'T') ); % get rid of all the unknowns
                I = unique(I);
                block_kmers_unique = block_kmers_unique(I,:);
                block_kmers_counts = uint32(block_kmers_counts(I));
                if(conservation_flag)
                    block_kmers_weights = block_kmers_weights(I);
                    block_kmers_weights_sqr = block_kmers_weights_sqr(I);
                end

                % take only kmers which have a certain prefix
                if(~isempty(prefix))
                    prefix_inds = strmatch(prefix, block_kmers_unique);
                    block_kmers_unique = block_kmers_unique(prefix_inds, :); block_kmers_counts = block_kmers_counts(prefix_inds);
                    if(conservation_flag)
                        block_kmers_weights = block_kmers_weights(prefix_inds);  block_kmers_weights_sqr = block_kmers_weights_sqr(prefix_inds);
                    end
                end
                doing_block = j;
                if(~isempty(block_kmers_unique)) % it could be that we didn't find any matching kmer
                    block_kmers_unique = pack_seqs(nt2int(block_kmers_unique)); % pack it and intersect with already established
                    if(conservation_flag)
                        [mid_block_kmers_unique mid_block_kmers_counts] = UnionWithCounts(mid_block_kmers_unique, mid_block_kmers_counts, ...
                            block_kmers_unique, [double(block_kmers_counts) block_kmers_weights block_kmers_weights_sqr]);
                    else
                        [mid_block_kmers_unique mid_block_kmers_counts] = UnionWithCounts(mid_block_kmers_unique, mid_block_kmers_counts, ...
                            block_kmers_unique, block_kmers_counts);
                    end
                    if (mod(j, mid_block_size) == 0) % || (j==num_blocks) )
                        doing_block = j;
                        [kmers_unique kmers_counts] = UnionWithCounts(kmers_unique, kmers_counts, ...
                            mid_block_kmers_unique, mid_block_kmers_counts);
                        mid_block_kmers_unique = []; mid_block_kmers_counts = [];
                    end
                end
            end % if exists file
        end
    end
    % unite again last time
    if(~isempty(mid_block_kmers_unique))
        [kmers_unique kmers_counts] = UnionWithCounts(kmers_unique, kmers_counts, ...
            mid_block_kmers_unique, mid_block_kmers_counts);
        mid_block_kmers_unique = []; mid_block_kmers_counts = [];
    end
    toc
    % save chromosome output - we don't want to save every block - just
    % at the end - maybe this should also move out of the loop on chromsomes
    if(save_blocks_flag == 0)
        if(conservation_flag) % save the data to output file for EACH chromosome, so we don't loose stuff. The problem is that currently this is done in a cumulative way
            if(~isempty(kmers_counts))
                kmers_weights = kmers_counts(:,2); kmers_weights_sqr = kmers_counts(:,3); kmers_counts = uint32(kmers_counts(:,1));
            else
                kmers_weights = []; kmers_weights_sqr = [];
            end
            save(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)], ...
                ['chr' chr_num2str(chr, organism_str) '_' output_file]), ...
                'kmers_unique', 'kmers_counts', 'kmers_weights', 'kmers_weights_sqr', 'num_blocks', 'regions_len', 'prefix');
        else
            save(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)], ...
                ['chr' chr_num2str(chr, organism_str) '_' output_file]), ...
                'kmers_unique', 'kmers_counts', 'num_blocks', 'regions_len', 'prefix');
        end
    end
end

% Now unite EVERYTHING at the end and save
tic;
if( (cumulative_flag == 0)&&(save_blocks_flag == 0) ) % in this case work on each chromosome seperately - faster and more parallelible but need to unite results at the end
    file_names = {};
    for i=1:length(chroms)
        chr = chroms(i);
        file_names{i} = fullfile(['chr' chr_num2str(chr, organism_str)], ['chr' chr_num2str(chr, organism_str) '_' output_file]);
    end
    Dummy = UniteKmersCountsFiles(fullfile(output_dir, ['k' num2str(k)]), file_names, ['all_' output_file], conservation_flag);
end
toc

