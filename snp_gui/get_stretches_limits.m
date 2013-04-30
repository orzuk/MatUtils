%function [samples_start_end_cell, samples_loh_cell, ret_snp_loc_vec, samples_big_aberr_vec] = get_stretches_limits(diag_table, samples)
function [samples_start_end_cell, samples_loh_cell, ret_snp_loc_vec, samples_big_aberr_vec] = get_stretches_limits(stretches_table, samples)

[C, IA, calls_columns] = intersect_order_by_first_gr(samples, stretches_table(1,:));
samples_start_end_cell = [];
samples_loh_cell = [];
ret_snp_loc_vec = [];
samples_big_aberr_vec = [];

snp_ind = strmatch('SNP_', stretches_table(:,1));
loc_column = strmatch(lower('LOCATION'), lower(stretches_table(1,:)));
chr_column = strmatch(lower('CHR'), lower(stretches_table(1,:)));
snp_chr_vec = cell_column_of_nums_to_mat(stretches_table(snp_ind, chr_column));
snp_loc_vec = cell_column_of_nums_to_mat(stretches_table(snp_ind, loc_column));

diag_ones_vec = zeros(1, max(snp_ind));
diag_ones_vec(snp_ind) = 1;
if(length(diag_ones_vec)>0)

    [diag_joined_ones_vec, diag_stretch_len_vec] = count_joined_ones(diag_ones_vec);
    num_stretches = length(diag_joined_ones_vec(1,:));

    stretch_ind_start_vec = zeros(1, num_stretches);
    stretch_ind_end_vec = zeros(1, num_stretches);
    for i = 1:num_stretches
        stretch_ind_start_vec(i) = diag_joined_ones_vec(2,i);
        stretch_ind_end_vec(i) = stretch_ind_start_vec(i)+diag_joined_ones_vec(1,i)-1;
    end
    % now count for each stretch
    num_samples = length(calls_columns);
    samples_start_end_cell = cell(num_samples, 4);
    samples_loh_cell = cell(num_samples, 1);
    samples_big_aberr_vec = zeros(num_samples, 1);
    loc_vec = [];
    temp = snp_ind;
    %     big_aberr_ind = [1:size(stretches_table, 1)];
    %     big_aberr_ind(temp) = [];
    for j = 1:num_stretches
        stretch_ind_start = stretch_ind_start_vec(j);
        stretch_ind_end = stretch_ind_end_vec(j);

        stertch_ind = [stretch_ind_start: stretch_ind_end];
        stretch_cell = stretches_table(stertch_ind, :);
        calls_cell = stretches_table(stertch_ind, calls_columns);
        stretch_loc_vec = cell_column_of_nums_to_mat(stretches_table(stertch_ind, loc_column));
        if(size(stretch_loc_vec,1) ~=1) stretch_loc_vec = stretch_loc_vec'; end
        loc_vec = [loc_vec stretch_loc_vec];
        for k = 1:num_samples
            sample_in_stretch_ind = strmatch('**', calls_cell(:,k));
            if(length(sample_in_stretch_ind)>0)
                temp_start = samples_start_end_cell{k,1};
                temp_end = samples_start_end_cell{k,2};
                stretch_len = samples_start_end_cell{k,3};
                [start_vec, end_vec] = ind_seq_into_start_end(sample_in_stretch_ind);
                temp_start = [temp_start stretch_loc_vec(start_vec)];
                temp_end = [temp_end stretch_loc_vec(end_vec)];
                samples_start_end_cell{k,1} = temp_start;
                samples_start_end_cell{k,2} = temp_end;
                % find stretch len
                for b = 1:length(start_vec)
                    ind = start_vec(b);
                    stretch_str = char(calls_cell(ind,k));
                    psik_ind = strfind(stretch_str,',');
                    stretch_len_str = stretch_str(psik_ind(1)+1:psik_ind(2)-1);
                    stretch_len = [stretch_len str2num(stretch_len_str)];
                end
                samples_start_end_cell{k,3} = stretch_len;
            end
            %now find loh/LOH
            temp = samples_loh_cell{k,1};
            loh_ind = strmatch('**LOH', calls_cell(:,k));
            if(length(loh_ind)>0)
                temp = [temp stretch_loc_vec(loh_ind)];
            end
            loh_ind = strmatch('loh', calls_cell(:,k));
            if(length(loh_ind)>0)
                temp = [temp stretch_loc_vec(loh_ind)];
            end
            samples_loh_cell{k,1} = temp;
        end
        %now find big aberrations - take the second stretch if there is (
        %before it there will always be big aberrations line

        if(num_stretches>1)
            stretch_ind_start = stretch_ind_start_vec(2);
        else
            stretch_ind_start = stretch_ind_start_vec(1);
        end
        if(stretch_ind_start>1)
            big_aberr_line = stretches_table(stretch_ind_start-1, calls_columns);
            amp_ind = strmatch('AMP', big_aberr_line);
            del_ind = strmatch('DEL', big_aberr_line);
            samples_big_aberr_vec(amp_ind)=1;
            samples_big_aberr_vec(del_ind)=-1;
        end
    end
    ret_snp_loc_vec = loc_vec;
end