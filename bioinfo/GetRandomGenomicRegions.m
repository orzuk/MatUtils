% Return coordinates to regions selected randomly and uniformly from the genome
% In the future we should enable randomization also from specific regions types, like promtoers, introns, exons etc.
% (For now this is still possible, but more cumbursome, and you have to specify this by requiring the input 
% regions to be only promoters/introns/etc.)
%
% Input:
% num_regions - # of regions to randomize
% regions_lengths - the lengths of the different regions
% genome_version - which genome version do we draw from (optional)
% input_chr_vec - sample only from these chroms (optional)
% input_pos_start_vec - sample only from these regions (optional)
% input_pos_end_vec - sample only from these regions (optional)
%
% Output:
% rand_regions - A struct with the random regions selected, with the fields
% chr_vec, pos_start_vec, pos_end_vec and strand. 
%
function rand_regions = GetRandomGenomicRegions(num_regions, regions_lengths, genome_version, ...
    input_chr_vec, input_pos_start_vec, input_pos_end_vec, varargin)

Assign24MammalsGlobalConstants; % get genome version, etc.

if(length(regions_lengths) == 1) % duplicate in case all regions are of the same length
    regions_lengths = regions_lengths + zeros(num_regions, 1);
end

if(~exist('input_pos_end_vec', 'var')) % here sample from the entire genome
    input_pos_end_vec = []; % make it empty
end
if(~isempty(input_pos_end_vec)) % here sample from given regions
    if(length(input_chr_vec) == 1) % duplicate scalar inputs
        input_chr_vec = repmat(input_chr_vec, num_regions, 1);
    end
    if(length(input_pos_start_vec) == 1) % duplicate scalar inputs
        input_pos_start_vec = repmat(input_pos_start_vec, num_regions, 1);
    end
    if(length(input_pos_end_vec) == 1) % duplicate scalar inputs
        input_pos_end_vec = repmat(input_pos_end_vec, num_regions, 1);
    end
    input_region_lengths = input_pos_end_vec - input_pos_start_vec; % compute lengths
    input_reg_inds = weighted_rand(vec2row(input_region_lengths), num_regions); % randomize regions
    rand_regions.chr_vec = input_chr_vec(input_reg_inds); % get chromosomes

    rand_regions.pos_start_vec = input_pos_start_vec(input_reg_inds) + ...
        ceil(rand(num_regions, 1) .* (input_region_lengths(input_reg_inds) - regions_lengths)); % randomize starts
    rand_regions.pos_end_vec = rand_regions.pos_start_vec + regions_lengths - 1; % get ends
else  % here sample from the entire genome
    [genome_len genome_chr_lens] = GetGenomeLength(genome_version); % get genome and chromosomes length
    cumsum_chr_lens = cumsum(genome_chr_lens);
    chr_starts_vec = [1 cumsum_chr_lens'+1];

    % We must not intersect the boundary between two chromosomes
    pos_start_vec = ceil(rand(num_regions, 1) .* genome_len);
    pos_end_vec = pos_start_vec + regions_lengths;  % compute the chroms ends

    bad_inds = 1;
    while(~isempty(bad_inds))
        bad_inds = [];
        for i=1:num_regions
            start_chr = find(pos_start_vec(i) >= chr_starts_vec, 1, 'last');
            end_chr = find(pos_end_vec(i) >= chr_starts_vec, 1, 'last');
            if(start_chr < end_chr)
                bad_inds = [bad_inds i];
            end
        end
        pos_start_vec(bad_inds) = ceil(rand(length(bad_inds), 1) .* genome_len);
        pos_end_vec(bad_inds) = pos_start_vec(bad_inds) + regions_lengths(bad_inds);  % compute the chroms ends
    end
    chroms = 1:organism_chr+2;
    rand_regions.chr_vec = zeros(num_regions, 1);
    rand_regions.pos_start_vec = zeros(num_regions, 1);
    rand_regions.pos_end_vec = zeros(num_regions, 1);

    for i=1:length(chroms)
        chr_inds = find( (pos_start_vec >= chr_starts_vec(i)) & ...
            (pos_start_vec < chr_starts_vec(i+1)) );
        rand_regions.chr_vec(chr_inds) = i;
        rand_regions.pos_start_vec(chr_inds) = pos_start_vec(chr_inds) - chr_starts_vec(i)+1;
        rand_regions.pos_end_vec(chr_inds) = pos_end_vec(chr_inds) - chr_starts_vec(i)+1;
    end
end

rand_regions.strand = round(rand(num_regions, 1)); % determine strand randomly


