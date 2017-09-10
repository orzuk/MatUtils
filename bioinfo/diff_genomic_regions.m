% Take all regions appearing in one set of genomic regions and not the other
%
% Input:
% chr_vec1 - chromosomes of first interval set (or file with first set of regions)
% pos_start_vec1 - starts of first interval set
% pos_end_vec1 - ends of first interval set
% chr_vec2 - chromosomes of second interval set (or file with second set of regions)
% pos_start_vec2 - starts of second interval set
% pos_end_vec2 - ends of second interval set
% genome_version - which genome we use (only relevant when using files)
% output_file - where to save output (optional). Can be .txt or .mat format
% diff_mode - do we take parts of intervals (if intervals intersect). Default is 'yes'
%
% Output:
% diff_chr_vec - chromosomes of intersection
% diff_start_vec - starts of intersection
% diff_end_vec - ends of intersection
% diff_inds1 - indices of intersect intervals in first intervals
%
% New: The input can be from a file. If you want the first regions to be
% from a file you should put a string with file name instead of chr_vec1,
% and then pos_start_vec1, pos_end_vec1 will not be used (so you can leave
% them empty). The same for chr_vec2 (so you can have one/two/zero of the
% intervals to come from a file). Note: if you read from a file you must
% specify the genome version (to enable correct conversion of X,Y chromosomes)
%
function [diff_chr_vec diff_start_vec diff_end_vec diff_inds1] = ...
    diff_genomic_regions(chr_vec1, pos_start_vec1, pos_end_vec1, ...
    chr_vec2, pos_start_vec2, pos_end_vec2, genome_version, output_file, diff_mode, varargin)

FULL_INTERVALS = 1; BREAK_INTERVALS = 0; % set difference mode 
if(~exist('diff_mode', 'var') || isempty(diff_mode))
    diff_mode = BREAK_INTERVALS; 
end

if(ischar(chr_vec1)) % first regions come from a file
%    file_name1 = chr_vec1;
    if(strcmp(suffix_from_file_name(chr_vec1), 'mat')) % .mat file
        load(chr_vec1); % assume format is: chr_vec, pos_start_vec, pos_end_vec
        chr_vec1 = chr_vec; pos_start_vec1 = pos_start_vec; pos_end_vec1 = pos_end_vec;
        clear chr_vec pos_start_vec pos_end_vec;
    else
        [chr_vec1 pos_start_vec1 pos_end_vec1] = ...
            ReadWigFileToMat(chr_vec1, genome_version); % here we must specify genome version
    end
else
%    file_name1 = '1st-regions';
end
if(ischar(chr_vec2)) % second regions come from a file
%    file_name2 = chr_vec2;
    if(strcmp(suffix_from_file_name(chr_vec2), 'mat')) % .mat file
        load(chr_vec2); % assume format is: chr_vec, pos_start_vec, pos_end_vec
        chr_vec2 = chr_vec; pos_start_vec2 = pos_start_vec; pos_end_vec2 = pos_end_vec;
        clear chr_vec pos_start_vec pos_end_vec;
    else
        [chr_vec2 pos_start_vec2 pos_end_vec2] = ...
            ReadWigFileToMat(chr_vec2, genome_version); % here we must specify genome version
    end
else
%    file_name2 = '2nd-regions';
end
chroms = intersect( unique(chr_vec1), unique(chr_vec2) );
diff_start_vec = []; diff_end_vec = []; diff_chr_vec = []; diff_inds1 = []; diff_inds2 = [];
for i=1:length(chroms)
%    chr_is = i
    chr_inds1 = find(chr_vec1 == chroms(i));
    chr_inds2 = find(chr_vec2 == chroms(i));
    [diff_start diff_end cur_diff_inds1 cur_diff_inds2] = ...
        intervals_diff(pos_start_vec1(chr_inds1), pos_end_vec1(chr_inds1), ...
        pos_start_vec2(chr_inds2), pos_end_vec2(chr_inds2), [], [], diff_mode);
    num_regions = length(diff_start);
    if(num_regions > 0)
        diff_start_vec = [diff_start_vec vec2row(diff_start)];
        diff_end_vec = [diff_end_vec vec2row(diff_end)];
        diff_chr_vec = [diff_chr_vec repmat(chroms(i), 1, num_regions)];
        diff_inds1 = [vec2row(diff_inds1) vec2row(chr_inds1(cur_diff_inds1))];
%        diff_inds2 = [vec2row(diff_inds2) vec2row(chr_inds2(diff_inds2))];
    end
end

if(exist('output_file', 'var') && (~isempty(output_file))) % save output
    if(strcmp(suffix_from_file_name(output_file), 'mat')) % save to .mat file
        save(output_file, 'diff_chr_vec', 'diff_start_vec', 'diff_end_vec', ...
            'diff_inds1', 'diff_inds2');
    else % save to .txt file
        save_regions_mat_as_text(output_file, diff_chr_vec, diff_start_vec, diff_end_vec, ...
            genome_version, [vec2column(diff_inds1) vec2column(diff_inds2)], 1);
    end
end


