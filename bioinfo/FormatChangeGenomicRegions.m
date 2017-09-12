% Convert .txt format chr5:start-end to numeric: 5 start end
% 
% Input: 
% regions_vec - vector with txt format
% genome_version (optional)
% 
% The output: 
% chr_vec - chromosome
% pos_start_vec - region starts
% pos_end_vec - region ends
%
function [chr_vec pos_start_vec pos_end_vec] = FormatChangeGenomicRegions(regions_vec, genome_version, varargin)

Assign24MammalsGlobalConstants;

n = size(regions_vec,1)
chr_vec = zeros(n,1); pos_start_vec = zeros(n,1); pos_end_vec = zeros(n,1);

for i=1:n
    if(mod(i,100) == 0)
        i_is = i
    end
    first_stop = strfind(regions_vec{i}, ':');
    second_stop = strfind(regions_vec{i}, '-');
    if(~isempty(strfind(regions_vec{i}, '_random'))) % get rid of random stuff
        chr_vec(i) = chr_str2num(regions_vec{i}(4:first_stop-8), organism_str );
    else
        chr_vec(i) = chr_str2num(regions_vec{i}(4:first_stop-1), organism_str );
    end
    pos_start_vec(i) = str2num(regions_vec{i}(first_stop+1:second_stop-1));
    pos_end_vec(i) = str2num(regions_vec{i}(second_stop+1:end));
end
