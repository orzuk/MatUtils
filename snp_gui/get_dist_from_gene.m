%function [dist_vec] = get_dist_from_gene(loc_start, loc_end, exon_start_vec, exon_end_vec, in_loc_vec)
function [dist_vec] = get_dist_from_gene(loc_start, loc_end, exon_start_vec, exon_end_vec, in_loc_vec)

num_loc = length(in_loc_vec);
dist_vec = zeros(num_loc, 1);
dist_vec(:) = -1;
dist_str_cell = cell(num_loc, 1);
dist_str_cell(:) = {''};
% find locations which are in/out gene
ret_vec = is_bet_vals_in_vec(loc_start, loc_end, in_loc_vec);
out_ind = find(ret_vec==0);
in_ind = find(ret_vec==1);
dist_vec(out_ind) = min(abs(in_loc_vec(out_ind)-loc_start), abs(in_loc_vec(out_ind)-loc_end));
% for i = 1:length(out_ind)
%     dist_str_cell{out_ind(i)} = num2str(dist_vec(out_ind(i)));
% end

% for each of in genes
num_in_genes = length(in_ind);
for i = 1:num_in_genes
    in_loc = in_loc_vec(in_ind(i));
    % is in utr
    min_exon_loc = min(min(exon_start_vec), min(exon_end_vec));
    max_exon_loc = max(max(exon_start_vec), max(exon_end_vec));
    if(is_bet_vals(min_exon_loc, max_exon_loc, in_loc))
        % exon or intron
        exon_ind = is_bet_vals(exon_start_vec, exon_end_vec, in_loc);
        exon_ind = find(exon_ind==1);
        if(length(exon_ind)>0)
            dist_str = 'CDS';
            dist_vec(in_ind(i)) = -3;
        else
            dist_str = 'intron';
            dist_vec(in_ind(i)) = -2;
        end
    else % in utr
        dist_str = 'UTR';
        dist_vec(in_ind(i)) = -1;
    end
%    dist_str_cell{in_ind(i)} = dist_str;
end

