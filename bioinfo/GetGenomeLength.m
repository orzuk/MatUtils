% Gets the length of all the chromosomes in a given genome. 
% for the first time we actually do the work, but later we just read the
% lengths from a file where they are stored. There is no input - the
% organism and genome version are defined in Assign24 ... 
% 
% Input: 
% genome_version - which genome version do we draw from (optional)
%
% Output: 
% genome_len - the length of the whole genome
% genome_chr_lens - a vector with the lengths of all chromosomes
%
function [genome_len genome_chr_lens] = GetGenomeLength(genome_version, varargin)

Assign24MammalsGlobalConstants;
data_dir = fullfile('../data', genome_version); 

num_chroms = organism_chr+2; % take also chr Y and Mitochondria 
if(~exist(fullfile(data_dir, 'genome_annotations.mat'), 'file'))
    genome_chr_lens = zeros(num_chroms,1);
    for i=1:num_chroms
        if(exist('Sequence', 'var'))
            clear Sequence 
        end
        do_chr_genome_len = i
        load(fullfile(data_dir, ['chr' chr_num2str(i, organism_str)]));
        genome_chr_lens(i) = length(Sequence);
    end
    save(fullfile(data_dir, 'genome_annotations.mat'), 'genome_chr_lens');
else
    load(fullfile(data_dir, 'genome_annotations.mat'));
end
genome_len = sum(genome_chr_lens);

