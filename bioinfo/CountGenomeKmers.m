% Go over the WHOLE genome and count the occurances of each kmer.
% This function handles the interfacing job but the counting itself is done
% by the function 'CountGenomicRegionsKmers' which is run through bsub from
% the farm, so we can divide the job and run in parralel
%
% Input:
% input_seq_dir - main directory with sequences in sub-directories
% input_scores_dir - directory with scores 
% output_kmers_dir - main output directory where to save kmers scores 
% organism_dir - sub-directory with organism genome version - sequence is here
% output_file - where to save results
% block_size - size of blocks used to divide work to blocks to avoid memory problems
% save_blocks_flag - 1: count kmers and generate blocks. 0: collect and combine blocks of kmers statistics
% chroms - which chromosomes to run
% k - kmer length
% pref_len - length of prefix by which to divide kmers to groups (4^pref_len groups overall)
% conservation_flag - new!! allow to keep conservation!!
%
% Output: saved to file output_file
%
function Dummy = CountGenomeKmers(input_seq_dir, input_scores_dir, output_kmers_dir, organism_dir, output_file, ...
    block_size, save_blocks_flag, cumulative_flag, chroms, k, pref_len, conservation_flag, clade_str, ...
    queue_str, mem_flag, in_matlab_flag) % last 3 parameters related to jobs scheduling

Assign24MammalsGlobalConstants;

if(~exist('conservation_flag', 'var') || isempty(conservation_flag))
    conservation_flag = 1; % temp
end
if(~exist('queue_str', 'var') || isempty(queue_str))
    queue_str = 'priority'; % temp
end
if(~exist('mem_flag', 'var') || isempty(mem_flag))
    mem_flag = 12; % temp: request more memory for uniting all blocks
end
if(~exist('in_matlab_flag', 'var') || isempty(in_matlab_flag))
    in_matlab_flag = 0;  % default: submit jobs
end
if(~exist('input_scores_dir', 'var'))
    input_scores_dir = []; 
end
if(~exist('output_kmers_dir', 'var') || isempty(output_kmers_dir))
    output_kmers_dir = input_seq_dir; 
end

if(save_blocks_flag)
    cumulative_flag = 0; % run chrom. by chrom.
else
    if(~exist('cumulative_flag', 'var') || isempty(cumulative_flag))
        cumulative_flag = 1; % unite all blocks to one file per chromosome 
    %    cumulative_flag = -1; % just unite stuff at the end
    end
end

if(pref_len > k) % this cannot be!!!
    changing_prefix_to = k
    pref_len = k;
end

output_dir = fullfile(output_kmers_dir, organism_dir);
prefix_kmers = all_kmers(pref_len); % generate all kmers of size pref_len (default we use is four)
num_prefix = 4^pref_len;

if(~exist(fullfile(output_dir, ['k' num2str(k)]), 'dir'))
    mkdir(fullfile(output_dir, ['k' num2str(k)]));
end
for i=1:length(chroms)
    chr = chroms(i);
    if(~exist( fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)]), 'dir' ))
        mkdir(fullfile(output_dir, ['k' num2str(k)], ['chr' chr_num2str(chr, organism_str)]));
    end
end
ctr=0;
for i=1:num_prefix
    if(pref_len == 0)
        cur_pref = '';
    else
        cur_pref = prefix_kmers(i,:);
    end
    prefix_output_file = [cur_pref, '_', output_file];
    looking_for_chr_file = fullfile(output_dir, ['k' num2str(k)], ['all_' prefix_output_file])
    if(~exist(fullfile(output_dir, ['k' num2str(k)], ['all_' prefix_output_file]), ...
            'file')) %%    if(~exist(prefix_output_file, 'file')) % only count if there isn't yet ..
        job_str = ['[A B C D E F] =  CountGenomicRegionsKmers([], [], [], ' num2str(k) ',' num2str(block_size) ',' num2str(save_blocks_flag) ...
            ',''' cur_pref ''', ' num2str(cumulative_flag) ', ' num2str(conservation_flag)  ', ''' ...
            clade_str  ''', ''' input_scores_dir ''', ''' input_seq_dir ''',''' organism_dir ''',''' output_dir ''',[' num2str(chroms) '], ''' prefix_output_file ''');']
        if(in_matlab_flag)
            eval(job_str);
        else
            SubmitMatlabJobToFarm(job_str, fullfile(output_dir, ['k' num2str(k)], ... % run job (.out file contains clade) 
                ['job_' prefix_output_file '_' clade_str '.out']), queue_str, [], [], mem_flag);
        end
        ctr=ctr+1
    end
end

Dummy = 0;


