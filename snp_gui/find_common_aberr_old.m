%function find_common_aberr(sample_names, user_dir, ChipName, snp_ids, chr_loc_vec, chr_vec, snp_gene_symbols, snps_gene_dist, del_amp_flag, thresh_del, thresh_amp, Q)
function [segments mat_genes_unique segments_samples sample_names chr_and_loc_and_num_samples] = ...
    find_common_aberr(sample_names, user_dir, handles, chip_snp_ids_ordered, chr_loc_vec, chr_vec, ...
    snp_gene_symbols, snp_gene_descr, del_amp_flag, Q, min_stretch, max_stretch ,ucsc_webpage, genes_db_struct, ...
    gene_symbols, gene_descr, hmm_out)


ChipName=handles.ChipName;
end_p_location=handles.end_p_location;

%%%%%%%%%%%%%%%%%
% for checkind dchip data
%[d_copy_num_mat, d_genotype_mat, d_snp_ids] = load_leukemia_dchip_copy_num(sample_names, ChipName);

[sample_pairs_cell, pairs_samples_ind] = samples_into_pairs(sample_names);


call_type_mat = zeros(size(hmm_out.genotype_mat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For samples with normal, move to call type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genotype_str_cell = cell(size(hmm_out.genotype_mat));
% put genotype paired types in paired samples
for i = 1:size(sample_pairs_cell, 1)
    %sample_ind = strmatch(char(sample_pairs_cell(i,2)), sample_names);
    call_type_vec = call_vecs_into_call_types(hmm_out.genotype_mat(:, pairs_samples_ind(i,1)), ...
        hmm_out.genotype_mat(:, pairs_samples_ind(i,2)));
    call_type_mat(:,pairs_samples_ind(i,2)) = call_type_vec;

    genotype_str_cell(:, pairs_samples_ind(i,2)) = calls_types_into_call_str(call_type_vec);
end
% put genotype types in non paired samples
samples_ind_no_pair = 1:length(sample_names);
if(length(sample_pairs_cell)>0)
    samples_ind_no_pair([pairs_samples_ind(:,1)' pairs_samples_ind(:,2)']) = [];
end
for i = 1:length(samples_ind_no_pair)
    genotype_str_cell(:, samples_ind_no_pair(i)) = calls_ind_into_call_str(...
        hmm_out.genotype_mat(:, samples_ind_no_pair(i)));
end
samples_ind_d = 1:length(sample_names);
if(length(sample_pairs_cell)>0)
    samples_ind_d(pairs_samples_ind(:,1)) = [];
end
%hmm_out.genotype_mat = hmm_out.genotype_mat(:, samples_ind_d);
hmm_out.average_copy_num_mat = hmm_out.average_copy_num_mat(:, samples_ind_d);
hmm_out.copy_num_mat = hmm_out.copy_num_mat(:, samples_ind_d);
call_type_mat = call_type_mat(:, samples_ind_d);
genotype_str_cell = genotype_str_cell(:, samples_ind_d);
sample_names = sample_names(samples_ind_d);
num_samples = length(sample_names);

% now create stretch mat
stretch_mat = zeros(size(hmm_out.average_copy_num_mat));
zero_one_mat = zeros(size(hmm_out.average_copy_num_mat));
if(del_amp_flag==1)% deletions
%    zero_one_mat(hmm_out.average_copy_num_mat<=thresh_del)=1;
    zero_one_mat = hmm_out.zero_one_mat_del(:,samples_ind_d);
elseif(del_amp_flag==2) % amplifications
%    zero_one_mat(hmm_out.average_copy_num_mat>=thresh_amp)=1;
    zero_one_mat = hmm_out.zero_one_mat_amp(:,samples_ind_d);
end

[Retention, Non_Inf, No_Call_both, No_Call_first, No_Call_second, LOH] = ...
    call_type_ind_into_call_type_str_ind();
if(del_amp_flag==1)% deletions
    zero_one_mat(call_type_mat == LOH) = 1;
    % now put zeros in retention (without a high ratio)
    zero_one_mat(call_type_mat == Retention) = 0;
end
num_chr = length(hmm_out.chr_num_snps_vec);

%prepare p and q arms for check_if_whole_chr_aberr
[dum idx]=ismember(hmm_out.data_snp_ids, chip_snp_ids_ordered);
hmm_out.SNP_loc=chr_loc_vec(idx);
snp_ind=1;
hmm_out.arm_num_snps_vec=zeros(1,2*num_chr);
for j = 1:num_chr
    chr_snp_ind = snp_ind:snp_ind+hmm_out.chr_num_snps_vec(j)-1;
    chr_snp_loc=hmm_out.SNP_loc(chr_snp_ind);
    hmm_out.arm_num_snps_vec(2*(j-1)+1)=length(find(chr_snp_loc<=end_p_location(j)));
    hmm_out.arm_num_snps_vec(2*(j-1)+2)=length(find(chr_snp_loc>end_p_location(j)));
    snp_ind = snp_ind + hmm_out.chr_num_snps_vec(j);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporary: flags for remove big aberrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remove_big_del_flag = 0;
remove_big_amp_flag = 1;

big_aberr_mat = zeros(num_samples, num_chr,2); % keeps information of which chr arms are amplified/deleted
for i = 1:num_samples
    sample = char(sample_names{i});
    snp_ind=1;
    big_aberr_mat_sample = zeros(1,num_chr*2);
    for j = 1:num_chr*2 %go over p and q arms
        if(hmm_out.arm_num_snps_vec(j)>0)
            arm_snp_ind = snp_ind:snp_ind+hmm_out.arm_num_snps_vec(j)-1;
            %remove whole-chromosome deletions\amplifications from the
            %analysis. libi: should take care of X chromosome: if it is
            %female should be 2, male - 1.
            [whole_chr_aberr_flag del_flag amp_flag] = check_if_whole_chr_aberr...
                (hmm_out.copy_num_mat(arm_snp_ind,i), del_amp_flag);
            if(del_flag)
                big_aberr_mat_sample(j) = -1;
            end
            if(amp_flag)
                big_aberr_mat_sample(j) = 1;
            end
            if(whole_chr_aberr_flag && (del_amp_flag==2 && remove_big_amp_flag) | ...
                    (del_amp_flag==1 && remove_big_del_flag) | (j== 23*2 | j==23*2-1))
%            if(whole_chr_aberr_flag)
                zero_one_mat(arm_snp_ind, i) = 0;
            else
                [joined_ones_vec, stretch_len_vec] = count_joined_ones(zero_one_mat(arm_snp_ind, i)');
                stretch_mat(arm_snp_ind, i) = stretch_len_vec';
            end
            snp_ind = snp_ind + hmm_out.arm_num_snps_vec(j);
        end
    end
    big_aberr_mat(i,:,:) = vec_into_mat(big_aberr_mat_sample,2)';
end
% Put large num in LOH so that it won't be removed when removing small
% stretches
strong_loh_flag = 0; % libi: when it is one the LOH get too much weight and the snps are scattered and not in stretches
stretch_mat_orig = stretch_mat;
if(del_amp_flag==1 && strong_loh_flag)
    for i = 1:num_samples
        loh_vec = zeros(size(call_type_mat(:,i)));
        loh_vec(call_type_mat(:,i) == LOH) = 1;
        stretch_mat(:, i) = snp_large_num_in_loh(loh_vec, stretch_mat(:, i));
    end
end

% libi: should insert this?
%    no_del_ind = calc_stretch_len_thresh_chr(del_mat(:,i), joined_ones_vec, frac_cut, ...
%         stretch_len_vec, loh_ind_mat(:, i), chr_num_vec);
%     stretch_mat(:,i) = stretch_len_vec';
%     stretch_mat_cut(:,i) = stretch_len_vec';




stretch_mat_cut = stretch_mat;
stretch_mat_cut(stretch_mat_cut>max_stretch) = max_stretch;
% remove small stretches, but leave those with LOH.
stretch_mat_cut(stretch_mat_cut<min_stretch) = 0;
stretch_mat_orig(stretch_mat_orig<min_stretch) = 0;

pvals_over = calc_integer_p_value_conv(stretch_mat_cut);
num_rejected = fdr(pvals_over, Q);
[p_val_sorted, ind_sorted] = sort(pvals_over);
if (isempty(num_rejected) || num_rejected==0 || isempty(num_rejected))
    num_rejected = 300;
end
stretch_snps = sort((ind_sorted(1:num_rejected)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare snps_data_cell for function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% libi: need to remove empty hmm_out.data_snp_ids
%num_snps = length(hmm_out.data_snp_ids);
[C, AI, BI] = intersect_order_by_first_gr(hmm_out.data_snp_ids, chip_snp_ids_ordered);
%if(length(C) == num_snps) libi: or should check why in hmm_out.data_snp_ids there
%are empty ids
%snp_id_cell = cell(num_snps, 5);

snps_data_cell(:, 1) = chip_snp_ids_ordered(BI);
snps_data_cell(:, 2) = num2cell(chr_vec(BI));
snps_data_cell(:, 3) = num2cell(chr_loc_vec(BI));
snps_data_cell(:, 4) = snp_gene_descr(BI);
snps_data_cell(:, 5) = snp_gene_symbols(BI);
if(del_amp_flag==1)
    del_amp_str = 'DEL';
end
if(del_amp_flag==2)
    del_amp_str = 'AMP';
end
allele_ratio_mat = zeros(size(hmm_out.average_copy_num_mat));
%snps_calls_cell = cell(size(hmm_out.average_copy_num_mat));
remove_big_aberr_flag = 1;
proj_name = 'SNP';
stretch_snp_cell = build_snp_stretch_annot_cell(...
    proj_name, del_amp_str, sample_names, remove_big_aberr_flag, ...
    stretch_mat_orig, stretch_snps, snps_data_cell, hmm_out.average_copy_num_mat, pvals_over, ...
    allele_ratio_mat, genotype_str_cell, big_aberr_mat, handles.end_p_location);

stretches_out_f = fullfile(user_dir,'output', [del_amp_str '_stretches_' ChipName '.xls']);
saveCellFile(stretch_snp_cell, stretches_out_f);


%%%%%%%%%% creating mat structure
snp_chr=stretch_snp_cell(2:end,2);
snp_loc=stretch_snp_cell(2:end,4);
snp_pval=stretch_snp_cell(2:end,5);
snp_samples=stretch_snp_cell(2:end,12);
empty_idx=find(cellfun('isempty',snp_loc));

%matrix with all SNPS
mat_all_snps=stretch_snp_cell(2:end,:);
mat_all_snps(empty_idx,:)=[];
genes_pval=cell2mat(mat_all_snps(:,5));
chr_and_loc_and_num_samples=cell2mat(mat_all_snps(:,[2 4 6])); %loc of SNP and num of samples del or amp
[genes_unique dum unique_idx]=unique(mat_all_snps(:,9));
[dum idx]=ismember(genes_unique, gene_symbols);
gene_descr_unique=gene_descr(idx);
mat_genes_unique=cell(length(genes_unique),12);

% 1: Gene
% 2:  Gene Descr.
% 3: Chr
% 4: Band
% 5: Location
% 6: P-value
% 7: Normal Aberr.
% 8: Genecards
% 9: UCSC
% 10: Stretches
% 11: Num samples involved
% 12: Stretches Length

for i=1:length(genes_unique)
    idx_gene = find(unique_idx==i);
    [dum idx_min_pval]=min(genes_pval(idx_gene));
    mat_genes_unique{i,1}=genes_unique{i};
    mat_genes_unique{i,2}=gene_descr_unique{i};
    mat_genes_unique(i,3:7)=mat_all_snps(idx_gene(idx_min_pval),[2:5 end]);
    mat_genes_unique{i,8}=['=HYPERLINK("http://www.genecards.org/cgi-bin/carddisp.pl?gene=' genes_unique{i} '", "Genecards")'];
    mat_genes_unique{i,9}=['=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg17&position=chr2&' ...
        'hgt.customText=' ucsc_webpage '/' [genes_unique{i} '_' ChipName '_ucsc.txt'] '", "Ucsc Link")'];
    idx=find(strcmp(stretch_snp_cell(:,9),genes_unique{i}));
    mat_genes_unique{i,10}=['=HYPERLINK("[' stretches_out_f ']a' num2str(idx(1)) '","Stretches")'];
    mat_genes_unique{i,11}=length(get_unique_list(mat_all_snps(idx_gene,12)));
    
    % stretches length
    idx=find(strcmp(genes_db_struct.gene_symbols,genes_unique{i}));
    gene_start=genes_db_struct.loc_start(idx(1));
    gene_end=genes_db_struct.loc_end(idx(1));
    chr=genes_db_struct.chr(idx(1));
    mat_genes_unique{i,12}=length(find(cell2mat(mat_all_snps(:,2))==chr & ...
        cell2mat(mat_all_snps(:,4)) >= gene_start & cell2mat(mat_all_snps(:,4)) <= gene_end));
end

empty_idx(end+1)=length(snp_loc)+1;
segments=zeros(length(empty_idx)-1,4);
segments_samples=cell(length(empty_idx)-1,1);
for i=1:length(empty_idx)-1
    segments(i,1)=snp_loc{empty_idx(i)+1};
    segments(i,2)=snp_loc{empty_idx(i+1)-1};
    if segments(i,1)==segments(i,2)
        segments(i,2)=segments(i,2)+1; %for display
    end
    segments(i,3)=log10(snp_pval{empty_idx(i)+1});
    segments(i,4)=snp_chr{empty_idx(i)+1};
    segments_samples{i}=get_unique_list(snp_samples(empty_idx(i)+1:empty_idx(i+1)-1));
end

if del_amp_flag==2
    segments(:,3)=-segments(:,3);
end



function samples_vec=get_unique_list(samples_vec)

if (length(samples_vec)==1 && isempty(samples_vec{1})) %|| isempty(samples_vec)
    samples_vec={};
    return;
end

samples_vec=strvcat(samples_vec);
samples_vec=samples_vec';
samples_vec=samples_vec(:);
samples_vec=textscan(samples_vec,'%s','delimiter',' ','multipledelimsasone',1);
samples_vec=unique(samples_vec{1});
