% Build snp annotations 
%function [stretch_snp_cell, stretch_vec_start, stretch_vec_end] = build_snp_stretch_annot_cell(proj_name, samples, remove_big_aberr_flag, stretch_snps_samples_mat, stretch_snps, stretch_snps_samples_mat, snps_data_cell, snps_hmm_cell, snps_calls_cell)
function [stretch_snp_cell, stretch_vec_start, stretch_vec_end] = build_snp_stretch_annot_cell...
    (proj_name, del_amp_str, samples, remove_big_aberr_flag, stretch_mat, stretch_snps,...
    snps_data_cell, raw_copy_num_mat, p_vals_vec, snps_calls_cell, ...
    big_aberr_mat, end_p_location)


load(fullfile('..','database','normal_variation.mat'));
if strcmp(del_amp_str,'DEL')
    freq_col=5;
else
    freq_col=6;
end

allele_ratio_flag = 0;
calls_flag = 1;
aberr_cell = {};
%aberr_cell = extract_big_aberr_file(proj_name);

if(size(samples,1)~=1) samples = samples'; end
res_cel_first_line = {'SNP ID' 'Chr' 'Band' 'Location' 'P-value' ['Num ' del_amp_str ' Samples'] ...
    'Num Samples involved (Stretch/loh/big aberr)' ' Num loh Samples' 'Assoc. Gene' ...
    'Gene Descr.' 'Link' 'Samples'};
num_columns_no_samples = length(res_cel_first_line);
res_cel_first_line_samples = samples;
res_cel_first_line = concat_cells(res_cel_first_line, res_cel_first_line_samples, 1);
res_cel_first_line=[res_cel_first_line [del_amp_str ' norm. freq.']]; %add normal variation data
gene_symbols_list = cell(0,1);
num_columns = length(res_cel_first_line);

stretch_mat_ones = stretch_mat;
stretch_mat_ones(find(stretch_mat>1)) = 1;

snp_stretch_sum = sum(stretch_mat_ones')';
num_stretch_snps = length(stretch_snps);

snp_id_cell = snps_data_cell(:, 1);
chr_vec = cell2mat(snps_data_cell(:, 2));
chr_loc_vec = cell2mat(snps_data_cell(:, 3));
genes_descr = snps_data_cell(:, 4);
gene_symbols = snps_data_cell(:, 5);
stretch_snp_cell = cell(num_stretch_snps, num_columns);
if(num_stretch_snps)
    stretch_snp_cell(:,1) = snp_id_cell(stretch_snps, 1);
    stretch_snp_cell(:,2) = num2cell(chr_vec(stretch_snps));
    stretch_snp_cell(:,3) = get_chr_locs_band_chrs(chr_vec(stretch_snps), chr_loc_vec(stretch_snps));
    stretch_snp_cell(:,4) = num2cell(chr_loc_vec(stretch_snps));
    stretch_snp_cell(:,5) = num2cell(p_vals_vec(stretch_snps));
    stretch_snp_cell(:,6) = num2cell(snp_stretch_sum(stretch_snps));
    stretch_snp_cell(:,9) = gene_symbols(stretch_snps, 1);
    stretch_snp_cell(:,10) = genes_descr(stretch_snps, 1);
    gene_symbols_list = concat_cells(gene_symbols_list, gene_symbols(stretch_snps, 1), 2);
    stretch_snp_cell(:,11) = create_gene_cards_link_column(gene_symbols(stretch_snps, 1));
    %put normal variation data
    
    for i=1:size(stretch_snp_cell,1)
        idx=find(stretch_snp_cell{i,4}>=normal_variation_data(:,3) & stretch_snp_cell{i,4}<=normal_variation_data(:,4));
        if ~isempty(idx)
            stretch_snp_cell{i,end}=normal_variation_data(idx(1),freq_col);
        end
    end
    
    % put stretch sample names
    sample_names_line_space = samples;
    space_line = cell(size(sample_names_line_space)); space_line(:) = {' '};
    sample_names_line_space = concat_cell_str_column(sample_names_line_space, space_line);
    samples_big_cell = repmat(sample_names_line_space, num_stretch_snps, 1);
    non_stretch_ind = find(stretch_mat_ones(stretch_snps,:) == 0);
    samples_big_cell(non_stretch_ind) = {''};
    stretch_samples_column = cell(num_stretch_snps, 1);
    stretch_samples_column(:) = {''};
    for s = 1:size(samples_big_cell, 2)
        stretch_samples_column = concat_cell_str_column(stretch_samples_column, samples_big_cell(:, s));
    end
    stretch_snp_cell(:,12) = stretch_samples_column;
    intensity_mat = zeros(num_stretch_snps, length(samples));
    num_samples = length(samples);
    for s = 1: num_samples% write intensity and calls for each sample
        psik_column = cell(num_stretch_snps, 1);
        psik_column(:) = {','};
        sample_stretch_ind = find(stretch_mat_ones(stretch_snps, s)>0);
        star_column = cell(length(sample_stretch_ind), 1);
        star_column(:) = {'**'};
        if(calls_flag)
            snps_call_types = snps_calls_cell(:,s);
            snps_call_types = lower(snps_call_types);
            snps_call_types(stretch_snps(sample_stretch_ind)) = concat_cell_str_column(star_column, ...
                upper(snps_call_types(stretch_snps(sample_stretch_ind))));
            call_intens_column = concat_cell_str_column(snps_call_types(stretch_snps,:), psik_column);
        else
            call_intens_column = cell(num_stretch_snps, 1);
            call_intens_column(:) = {''};
        end
        raw_copy_vec = raw_copy_num_mat(:,s);
        intensity_mat(:,s) = raw_copy_vec(stretch_snps);
        % add stretch length of each sample
        if(~calls_flag)
            call_intens_column(sample_stretch_ind) = concat_cell_str_column(star_column, ...
                call_intens_column(sample_stretch_ind));
        end
        call_intens_column = concat_cell_str_column(call_intens_column, ...
            cellstr(num2str(stretch_mat(stretch_snps,s))));
        call_intens_column = concat_cell_str_column(call_intens_column, psik_column);
        call_intens_column = concat_cell_str_column(call_intens_column, ...
            cellstr(num2str(raw_copy_vec(stretch_snps,:))));
        if(allele_ratio_flag)
            call_intens_column = concat_cell_str_column(call_intens_column, psik_column);
            call_intens_column = concat_cell_str_column(call_intens_column, ...
                cellstr(num2str(allele_ratio_mat(stretch_snps,s))));
        end
        stretch_snp_cell(:,num_columns_no_samples+s) = call_intens_column;
    end
    %now look on stretches
    dist_next_stretch = stretch_snps(2:end)-stretch_snps(1:end-1);
    % put a sign if the chr changes
    chr_num_diff = chr_vec(stretch_snps(2:end))-chr_vec(stretch_snps(1:end-1));
    if(size(dist_next_stretch,1) ~= size(chr_num_diff, 1)) dist_next_stretch = dist_next_stretch'; end
    dist_next_stretch(find(chr_num_diff>0 & dist_next_stretch==1)) = 2;
    jump_ind = find(dist_next_stretch>1);
    num_jump_ind = length(jump_ind);
    % first stretch
    first_ind = stretch_snps(1);
    if(num_jump_ind)
        end_ind = stretch_snps(jump_ind(1));
    else
        end_ind = stretch_snps(end);
    end
    %     stretch_vec_start = first_ind;
    %     stretch_vec_end = end_ind;
    stretch_vec_start = [];
    stretch_vec_end = [];
    % put loh number
    loh_ind_mat = zeros(size(stretch_snp_cell,1), num_samples);
    loh_ind = strmatch('loh', stretch_snp_cell(:, num_columns_no_samples+1:num_columns_no_samples+num_samples));
    loh_ind_mat(loh_ind) = 1;
    sum_loh = sum(loh_ind_mat')';
    stretch_snp_cell(:,8) = num2cell(sum_loh);
    big_aberr_ind_mat = zeros(size(stretch_snp_cell,1), num_samples); % fill it in the next loop
    big_aberr_ind_mat_curr_stretch_start = 1;
    big_aberr_ind_mat_curr_stretch_end = 1;
    if(num_jump_ind)
        for b = 0:num_jump_ind
            if(b==0)
                curr_jump_ind = 1;
                curr_jump_ind_plus1 = curr_jump_ind;
            else
                curr_jump_ind = jump_ind(b);% the jump is from curr_jump_ind to curr_jump_ind+1
                curr_jump_ind_plus1 = curr_jump_ind+1;
            end
            dist_stretch_snp_cell = cell(1, num_columns); dist_stretch_snp_cell(:) = {''};
            num_gap_snps = dist_next_stretch(curr_jump_ind)-1;
            dist_stretch_snp_cell{1,1} = ['- ' num2str(num_gap_snps) ' -'];
            if(b < num_jump_ind)
                next_jump_ind = jump_ind(b+1);
            else
                next_jump_ind = num_stretch_snps;
            end

            stretch_vec_start = [stretch_vec_start stretch_snps(curr_jump_ind_plus1)];
            stretch_vec_end = [stretch_vec_end stretch_snps(next_jump_ind)];
            stretch_ind = find(stretch_mat_ones(stretch_snps(curr_jump_ind_plus1:next_jump_ind),:) >0);
            curr_intensity_mat = intensity_mat(curr_jump_ind_plus1:next_jump_ind,:);
            stretch_len = size(curr_intensity_mat, 1);
            big_aberr_ind_mat_curr_stretch_end = big_aberr_ind_mat_curr_stretch_start+stretch_len-1;
            mean_itens = mean(curr_intensity_mat(stretch_ind));
            %            if(remove_big_aberr_flag)
            curr_chr_loc_vec = chr_loc_vec(stretch_snps(curr_jump_ind_plus1:next_jump_ind));
            chr_num = chr_vec(stretch_snps(curr_jump_ind_plus1));
            for s = 1:length(samples)
%                check if this chromosome is in big aberrations in this sample
                big_aberr_mat_sample = zeros(24,2);
                big_aberr_mat_sample(:) = big_aberr_mat(s,:,:);
                if(chr_num==23)
                    t=1;
                end
                [ret, aberr_str] = check_if_chr_loc_big_aberr(big_aberr_mat_sample, end_p_location, chr_num, curr_chr_loc_vec);
                if(ret)
                    dist_stretch_snp_cell{1,num_columns_no_samples+s} = aberr_str;
                end
            end
            %           end
            %            dist_stretch_snp_cell{1,2} = [num2str(mean_itens) ':Avg. Raw Copy num' ];

            if(b==0)
                stretch_snp_cell1 = stretch_snp_cell([], :);
                stretch_snp_cell2 = stretch_snp_cell(1:end, :);
            else %libi: check this,
                %curr_jump_ind_in_cell = curr_jump_ind+b-1;
                curr_jump_ind_in_cell = curr_jump_ind+b;
                stretch_snp_cell1 = stretch_snp_cell(1:curr_jump_ind_in_cell, :);
                stretch_snp_cell2 = stretch_snp_cell(curr_jump_ind_in_cell+1:end, :);
            end
            % count loh samples, big aberr samples
            samples_ind = num_columns_no_samples+1:num_columns_no_samples+num_samples;
            big_aberr_line = dist_stretch_snp_cell(1,samples_ind);
            num_samples = length(samples);
            %stretch_samples_big_aberr_ind = get_big_aberr_samples_ind(big_aberr_line, del_amp_str);
            stretch_samples_big_aberr_ind = [];
            stretch_snp_cell = concat_cells(stretch_snp_cell1, dist_stretch_snp_cell, 2);
            stretch_snp_cell = concat_cells(stretch_snp_cell, stretch_snp_cell2, 2);
            big_aberr_ind_mat(big_aberr_ind_mat_curr_stretch_start:big_aberr_ind_mat_curr_stretch_end,...
                stretch_samples_big_aberr_ind) = 1;
            big_aberr_ind_mat_curr_stretch_start = big_aberr_ind_mat_curr_stretch_end+1;

        end
    end
    % now add the number of samples involved.
    cell_snp_ind = strmatch('SNP_A', stretch_snp_cell(:,1));
    inv_samples = big_aberr_ind_mat | loh_ind_mat;
    num_inv_samples = sum(inv_samples')';
    num_del_samples = cell2mat(stretch_snp_cell(cell_snp_ind,6));
    stretch_snp_cell(cell_snp_ind,7) = num2cell(num_inv_samples+num_del_samples);
    num_big_aberr_samples = sum(big_aberr_ind_mat')';
    %stretch_snp_cell(cell_snp_ind,9) = num2cell(num_big_aberr_samples);
    stretch_snp_cell = concat_cells(res_cel_first_line, stretch_snp_cell, 2);
end


