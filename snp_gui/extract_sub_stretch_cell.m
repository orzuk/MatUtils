%function [sub_stretch_cell, samples] =
%extract_sub_stretch_cell(stretch_cell, chr_num, loc_start, loc_end)
function [sub_stretch_cell, samples] = extract_sub_stretch_cell(stretch_cell, chr_num, loc_start, loc_end)

sub_stretch_cell = [];
samples = [];
if(length(stretch_cell))
    first_line = stretch_cell(1,:);
    sub_stretch_cell = first_line;
    chr_column = min(strmatch('Chr', first_line));
    loc_column = min(strmatch('Location', first_line));
    samples_column = min(strmatch('Samples', first_line));
    samples = first_line(samples_column+1:end-1)';

    snp_ind = strmatch('SNP_', stretch_cell(:,1));

    if(length(snp_ind))
        chr_vec = cell_column_of_nums_to_mat(stretch_cell(snp_ind, chr_column));
        loc_vec = cell_column_of_nums_to_mat(stretch_cell(snp_ind, loc_column));

        chr_ind = find(chr_vec==chr_num);
        if(length(chr_ind)>0)
            min_loc = min(find(loc_vec(chr_ind) >= loc_start));
            max_loc = max(find(loc_vec(chr_ind) <= loc_end));
            if(length(min_loc)>0 & length(max_loc)>0)

                start_ind = snp_ind(chr_ind(min_loc));
                % find the last empty line before start_ind - for big aberrations
                % info
                empty_inds = strmatch('-', stretch_cell(:,1));
                last_empty_ind = max(empty_inds(find(empty_inds<start_ind)));
                end_ind = snp_ind(chr_ind(max_loc));
                sub_stretch_cell2 = stretch_cell(start_ind:end_ind, :);
                if(length(last_empty_ind)>0)
                    sub_stretch_cell2 = concat_cells(stretch_cell(last_empty_ind,:), sub_stretch_cell2, 2);
                end
                sub_stretch_cell = concat_cells(sub_stretch_cell, sub_stretch_cell2, 2);
            end
        end
    end
end