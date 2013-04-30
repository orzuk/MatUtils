% Intersect two sets of genomic regions
%
% Input:
% chr_vec1 - chromosomes of first interval set (or file with first set of regions)
% pos_start_vec1 - starts of first interval set
% pos_end_vec1 - ends of first interval set
% chr_vec2 - chromosomes of second interval set (or file with second set of regions)
% pos_start_vec2 - starts of second interval set
% pos_end_vec2 - ends of second interval set
% genome_version - which genome we use (only relevant when using files)
% venn_plot_flag - plot a venn diagram of intersection
% output_file - where to save output (optional). Can be .txt or .mat format
% fig_outfile - where to save the figure (if plotting a venn diagram)
% venn_str - ?? 
% 
% Output:
% intersect_chr_vec - chromosomes of intersection
% intersect_start_vec - starts of intersection
% intersect_end_vec - ends of intersection
% intersect_inds1 - indices of intersect intervals in first intervals
% intersect_inds2 - indices of intersect intervals in second intervals
%
% New: The input can be from a file. If you want the first regions to be
% from a file you should put a string with file name instead of chr_vec1,
% and then pos_start_vec1, pos_end_vec1 will not be used (so you can leave
% them empty). The same for chr_vec2 (so you can have one/two/zero of the
% intervals to come from a file). Note: if you read from a file you must
% specify the genome version (to enable correct conversion of X,Y chromosomes)
%
function [intersect_chr_vec intersect_start_vec intersect_end_vec intersect_inds1 intersect_inds2 ...
    regions_len1 regions_len2 regions_inter_len] = ...
    intersect_genomic_regions(chr_vec1, pos_start_vec1, pos_end_vec1, ...
    chr_vec2, pos_start_vec2, pos_end_vec2, genome_version, venn_plot_flag, ...
    output_file, fig_outfile, venn_str, varargin)

if(~exist('fig_outfile', 'var'))
    fig_outfile = [];
end

if(ischar(chr_vec1)) % first regions come from a file
    file_name1 = chr_vec1;
    if(strcmp(suffix_from_file_name(chr_vec1), 'mat')) % .mat file
        load(chr_vec1, 'chr_vec', 'pos_start_vec', 'pos_end_vec'); % assume format is: chr_vec, pos_start_vec, pos_end_vec
        chr_vec1 = chr_vec; pos_start_vec1 = pos_start_vec; pos_end_vec1 = pos_end_vec;
        clear chr_vec pos_start_vec pos_end_vec;
    else
        [chr_vec1 pos_start_vec1 pos_end_vec1] = ...
            ReadWigFileToMat(chr_vec1, genome_version); % here we must specify genome version
    end
else
    file_name1 = '1st-regions';
end
if(ischar(chr_vec2)) % second regions come from a file
    file_name2 = chr_vec2;
    if(strcmp(suffix_from_file_name(chr_vec2), 'mat')) % .mat file
        load(chr_vec2, 'chr_vec', 'pos_start_vec', 'pos_end_vec'); % assume format is: chr_vec, pos_start_vec, pos_end_vec
        chr_vec2 = chr_vec; pos_start_vec2 = pos_start_vec; pos_end_vec2 = pos_end_vec;
        clear chr_vec pos_start_vec pos_end_vec;
    else
        [chr_vec2 pos_start_vec2 pos_end_vec2] = ...
            ReadWigFileToMat(chr_vec2, genome_version); % here we must specify genome version
    end
else
    file_name2 = '2nd-regions';
end
chroms = intersect( unique(chr_vec1), unique(chr_vec2) );
intersect_start_vec = []; intersect_end_vec = []; intersect_chr_vec = []; intersect_inds1 = []; intersect_inds2 = [];
for i=1:length(chroms)
    chr_inds1 = find(chr_vec1 == chroms(i));
    chr_inds2 = find(chr_vec2 == chroms(i));
    [inter_start inter_end inter_inds1 inter_inds2] = ...
        intervals_intersect(pos_start_vec1(chr_inds1), pos_end_vec1(chr_inds1), ...
        pos_start_vec2(chr_inds2), pos_end_vec2(chr_inds2));
    num_regions = length(inter_start);
    if(num_regions > 0)
        intersect_start_vec = [intersect_start_vec vec2row(inter_start)];
        intersect_end_vec = [intersect_end_vec vec2row(inter_end)];
        intersect_chr_vec = [intersect_chr_vec repmat(chroms(i), 1, num_regions)];
        intersect_inds1 = [intersect_inds1 vec2row(chr_inds1(inter_inds1))];
        intersect_inds2 = [intersect_inds2 vec2row(chr_inds2(inter_inds2))];
    end
end

regions_len1 = []; regions_len2 = []; regions_inter_len = []; 
if(exist('venn_plot_flag', 'var') && ~isempty(venn_plot_flag))
    if(venn_plot_flag)
        [regions_len1 regions_len2 regions_inter_len] = ...
            venn_genomic_regions(chr_vec1, pos_start_vec1, pos_end_vec1, remove_dir_from_file_name(file_name1), ...
            chr_vec2, pos_start_vec2, pos_end_vec2, remove_dir_from_file_name(file_name2), ...
            intersect_chr_vec, intersect_start_vec, intersect_end_vec, 'intersection', genome_version, ...
            1, fig_outfile)
    end
end
if(exist('output_file', 'var') && (~isempty(output_file))) % save output
    if(strcmp(suffix_from_file_name(output_file), 'mat')) % save to .mat file
        save(output_file, 'intersect_chr_vec', 'intersect_start_vec', 'intersect_end_vec', ...
            'intersect_inds1', 'intersect_inds2');
    else % save to .txt file
        save_regions_mat_as_text(output_file, intersect_chr_vec, intersect_start_vec, intersect_end_vec, ...
            genome_version, [vec2column(intersect_inds1) vec2column(intersect_inds2)], 1);
    end
end


