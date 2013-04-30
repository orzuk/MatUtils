% Stratify a set of genomic positions into different genomic regions.
%
% Input:
% chr_vec - the chromosomes
% pos_start_vec - the positions (starts or all)
% pos_end_vec - positions ends (optional)
% genome_version - the version of genome used 
% genomic_regions - a vector of starts and ends of regions for each chromosome
% is_sorted - flag saying if we need to sort the regions. (optional - default is 'yes')
%
% Output:
% reg_types_vec - for each position its type is returned
% reg_chr_vec - chroms of the returned type elements 
% reg_start_vec - starts of the returned type elements
% reg_end_vec - ends of the returned type elements
% reg_inds_vec - indices of the returned type elements in the input elements
% reg_types_counts_vec - histogram of each type in the regions (normalized)
%
function  [reg_types_vec reg_chr_vec reg_pos_start_vec reg_pos_end_vec reg_inds_vec reg_types_counts_vec] = ...
    StratifyToGenomicRegions(chr_vec, pos_start_vec, pos_end_vec, ...
     genome_version, genomic_regions, is_sorted, varargin)

if(~exist('pos_end_vec', 'var') || (isempty(pos_end_vec) && ~isempty(pos_start_vec)))  % Old function: work on single positions
    reg_types_vec = zeros(1, length(pos_start_vec), 'single')-1; % -1 indicates that we don't know yet the region type
    chroms = unique(chr_vec);
    if(size(chr_vec) == 1) % in this case there is only one chromosome and we indicate it by a number
        chr_vec = repmat(chr_vec, length(pos_start_vec), 1);
    end
    if(~exist('is_sorted', 'var') || isempty(is_sorted))
        is_sorted = 1;
    end
    if(~is_sorted)
        size_chr = size(chr_vec)
        size_pos = size(pos_start_vec)
        [pos_start_vec sort_perm] = sort(pos_start_vec); % we don't assume it's sorted
        chr_vec = chr_vec(sort_perm);
    end

    for i=1:length(chroms)
        chr = chroms(i);  stratify_chr = chr
        chr_inds = find(chr_vec == chr);
        chr_pos_start_vec = pos_start_vec(chr_inds);  % no need to sort positions (assume regions are already sorted)
        start_pos = min(chr_pos_start_vec); end_pos = max(chr_pos_start_vec);

        if( (~isempty(genomic_regions.chr_pos_start_vec{chr})) && (~isempty(start_pos)) )
            inter_inds = intersect( find(genomic_regions.chr_pos_start_vec{chr} <= end_pos), ...
                find(genomic_regions.chr_pos_end_vec{chr} >= start_pos) ); % look only at relevant indices
        else
            inter_inds = [];
        end
        for j=inter_inds'    %       1:  length(genomic_regions.chr_pos_start_vec{chr}) % loop over all segments
            I = find(chr_pos_start_vec >= genomic_regions.chr_pos_start_vec{chr}(j), 1);
            J = find(chr_pos_start_vec <= genomic_regions.chr_pos_end_vec{chr}(j), 1, 'last');
            reg_types_vec(chr_inds(I:J)) = genomic_regions.chr_region_type_vec{chr}(j);  % update region type for intersection
        end
        %    for j=1:length(genomic_regions.chr_pos_start_vec{chr}) % loop over all segments
        %        [cur_region_pos I J]= intersect(chr_pos_start_vec, ...
        %            genomic_regions.chr_pos_start_vec{chr}(j):genomic_regions.chr_pos_end_vec{chr}(j)); % perform (costly inefficient) intersection
        %        reg_types_vec(chr_inds(I)) = genomic_regions.chr_region_type_vec{chr}(j);  % update region type for intersection
        %    end
    end

    if(~is_sorted) % sort back
        reg_types_vec = reg_types_vec(inv_perm(sort_perm));
    end
else  % new function: work on intervals
    Assign24MammalsGlobalConstants;

    if(~exist('genomic_regions', 'var') || isempty(genomic_regions))
        genomic_regions = load(genome_annotations_file); % load annotation to perform intersect
    end
    chroms = unique(chr_vec);
    if(size(chr_vec) == 1) % in this case there is only one chromosome and we indicate it by a number
        chr_vec = repmat(chr_vec, length(pos_start_vec), 1);
    end
    reg_chr_vec = []; reg_pos_start_vec = []; reg_pos_end_vec = []; reg_types_vec = []; reg_inds_vec = [];
    for i=1:length(chroms) % loop on chromes
        chr = chroms(i);
        chr_inds = find(chr_vec == chr);
        all_regions_types = unique(genomic_regions.chr_region_type_vec{chr});
        for r=1:length(all_regions_types)
            chr_reg_inds = find(genomic_regions.chr_region_type_vec{chr} == all_regions_types(r));
            [chr_reg_intervals_start_pos chr_reg_intervals_end_pos chr_reg_intervals_inds1 chr_reg_intervals_inds2] = ...
                intervals_intersect(pos_start_vec(chr_inds), pos_end_vec(chr_inds), ...
                genomic_regions.chr_pos_start_vec{chr}(chr_reg_inds), genomic_regions.chr_pos_end_vec{chr}(chr_reg_inds));
            num_regs = length(chr_reg_intervals_start_pos);

            reg_chr_vec = [reg_chr_vec repmat(chr, 1, num_regs)];
            reg_pos_start_vec = [reg_pos_start_vec vec2row(chr_reg_intervals_start_pos)];
            reg_pos_end_vec = [reg_pos_end_vec vec2row(chr_reg_intervals_end_pos)];
            reg_types_vec = [reg_types_vec repmat(all_regions_types(r), 1, num_regs)];
            reg_inds_vec = [reg_inds_vec vec2row(chr_reg_intervals_inds1)];
        end
    end
end

reg_types_counts_vec = zeros(length(genomic_types_vec),1);
for i = 1:length(genomic_types_vec)
    cur_type_inds = find(reg_types_vec == genomic_types_vec(i));
    reg_types_counts_vec(i) = sum( reg_pos_end_vec(cur_type_inds) - reg_pos_start_vec(cur_type_inds) + 1);
end
reg_types_counts_vec = normalize(reg_types_counts_vec); % normalize



