% Merge (and remove overlap) a set of genomic regions
%
% Input:
% chr_vec - chromosomes of intervals
% pos_start_vec - intervals starting positions
% pos_end_vec - intervals ending positions
% genome_version - which genome we use (only relevant when using files)
% output_file - where to save output (optional). Can be .txt or .mat format
%
% Output:
% merge_chr_vec - chromosomes of merged regions
% merge_start_vec - merged starting positions
% merged_end_vec - merged ending positions
% merge_inds_vec - indices of original intervals in merged intervals
%
% New: The input can be from a file. If you want the regions to be
% from a file you should put a string with file name instead of chr_vec,
% and then pos_start_vec, pos_end_vec will not be used (so you can leave
% them empty). Note: if you read from a file you must
% specify the genome version (to enable correct conversion of X,Y chromosomes)
%
function [merge_chr_vec merge_start_vec merge_end_vec merge_inds_vec] = ...
    merge_genomic_regions(chr_vec, pos_start_vec, pos_end_vec, genome_version, output_file, varargin)

% inside_merging = 1111
% chr_vec_is = chr_vec
if(iscell(chr_vec)) % New! enable a list of files - need to read them one by one
%     inside_if_cell = 99999
    file_names = chr_vec;
    chr_vec = []; pos_start_vec = []; pos_end_vec = [];
    for i=1:length(file_names)
        if(exist(file_names{i}, 'file'))  % load only existing files
%             load_file = i
%            file_is = file_names{i}
            cur_regions = load(file_names{i});
            chr_vec = [chr_vec vec2row(cur_regions.chr_vec)];
            pos_start_vec = [pos_start_vec vec2row(cur_regions.pos_start_vec)];
            pos_end_vec = [pos_end_vec vec2row(cur_regions.pos_end_vec)];
            cur_len = length(chr_vec);
        end
    end
end
if(ischar(chr_vec)) % first regions come from a file
    if(strcmp(suffix_from_file_name(chr_vec1), 'mat')) % .mat file
        load(chr_vec); % assume format is: chr_vec, pos_start_vec, pos_end_vec
    else
        [chr_vec pos_start_vec pos_end_vec] = ...
            ReadWigFileToMat(chr_vec, genome_version); % here we must specify genome version
    end
end
chroms = unique(chr_vec);
merge_chr_vec = []; merge_start_vec = []; merge_end_vec = []; merge_inds_vec = [];
for i=1:length(chroms)
    chr_inds = find(chr_vec == chroms(i));
    merge_chr = chroms(i)
    merge_len = length(chr_inds)
    [merge_start merge_end merge_inds] = ...
        intervals_merge(pos_start_vec(chr_inds), pos_end_vec(chr_inds));
    num_regions = length(merge_start);
    if(num_regions > 0)
        merge_start_vec = [merge_start_vec vec2row(merge_start)];
        merge_end_vec = [merge_end_vec  vec2row(merge_end)];
        merge_chr_vec = [merge_chr_vec repmat(chroms(i), 1, num_regions)];
        merge_inds_vec = [merge_inds_vec  vec2row(chr_inds(merge_inds))];
    end
end
chr_vec = merge_chr_vec; % change naming back 
pos_start_vec = merge_start_vec; 
pos_end_vec = merge_end_vec;
inds_vec = merge_inds_vec;
if(exist('output_file', 'var')) % save output
    if(strcmp(suffix_from_file_name(output_file), 'mat')) % save to .mat file
        save(output_file, 'chr_vec', 'pos_start_vec', 'pos_end_vec', 'inds_vec');
    else % save to .txt file
        save_regions_mat_as_text(output_file, chr_vec, pos_start_vec, pos_end_vec, ...
            genome_version, vec2column(inds_vec), 1);
    end
end



