%function zero_one_mat = put_one_in_putative_aberr(del_amp_flag, raw_copy_num_mat, chr_num_snps_vec, copy_num_mat, thresh_del, thresh_amp, gender, chr_x_ind)
function zero_one_mat = put_one_in_putative_aberr(del_amp_flag, raw_copy_num_mat, chr_num_snps_vec, copy_num_mat, thresh_del, thresh_amp, gender, chr_x_ind)
 
% amp_quantile = 0.85;
% del_quantile = 0.35;
% amp_quantile = 0.9;
% del_quantile = 0.15;
% amp_quantile = 0.85;
% del_quantile = 0.25;
amp_quantile = 0.8;
del_quantile = 0.35;

zero_one_mat = zeros(size(copy_num_mat), 'single');
% if(del_amp_flag==1)% deletions
%     zero_one_mat(average_copy_num_mat<=thresh_del)=1;
% elseif(del_amp_flag==2) % amplifications
%     zero_one_mat(average_copy_num_mat>=thresh_amp)=1;
% end
two_copy_on_copy_diff = 0.8;
smooth_flag = 1;
%smooth_flag = 0;

num_snps = size(zero_one_mat,1);
non_chr_x_ind = [1:num_snps];
non_chr_x_ind(chr_x_ind) = [];
num_samples = size(copy_num_mat, 2);
num_chr = min(length(chr_num_snps_vec), 22); % don't run on X chromosome
for i = 1:num_samples % for each sample find copy number 2 snps and according to them decide on threshold
    diploid_snps_vec = zeros(num_snps,1);
    snp_ind = 1;
    raw_copy_num_vec_smooth = raw_copy_num_mat(:,i);
    for j = 1:num_chr
        if(chr_num_snps_vec(j)>0)
            chr_snp_ind = [snp_ind:snp_ind+chr_num_snps_vec(j)-1];
            if(smooth_flag)
                raw_copy_num_vec_smooth(chr_snp_ind) = smooth(raw_copy_num_vec_smooth(chr_snp_ind));
            end
            dip_snps_chr = find(copy_num_mat(chr_snp_ind,i)==2);
            frac_diploid = length(dip_snps_chr)/length(chr_snp_ind);
            if(frac_diploid >0.7)
                diploid_snps_vec(chr_snp_ind(dip_snps_chr)) = 1;
            end
            snp_ind = snp_ind + chr_num_snps_vec(j);
        end
    end
    if(del_amp_flag==1)% deletions
%        s_del_thresh = quantile(raw_copy_num_mat(find(diploid_snps_vec==1),i), del_quantile);
        s_del_thresh = quantile(raw_copy_num_vec_smooth(find(diploid_snps_vec==1)), del_quantile);
        if(isnan(s_del_thresh))
            s_del_thresh = thresh_del;
        end
        if (strcmpi(char(gender{i}),'M'))
%             zero_one_mat(chr_x_ind(raw_copy_num_mat(chr_x_ind,i)<s_del_thresh-two_copy_on_copy_diff),i) = 1;
%             zero_one_mat(non_chr_x_ind(raw_copy_num_mat(non_chr_x_ind,i)<s_del_thresh),i) = 1;
            zero_one_mat(chr_x_ind(raw_copy_num_vec_smooth(chr_x_ind)<s_del_thresh-two_copy_on_copy_diff),i) = 1;
            zero_one_mat(non_chr_x_ind(raw_copy_num_vec_smooth(non_chr_x_ind)<s_del_thresh),i) = 1;
        else
%            zero_one_mat(raw_copy_num_mat(:,i)<s_del_thresh,i) = 1;
            zero_one_mat(raw_copy_num_vec_smooth<s_del_thresh, i) = 1;
        end
    else % amplifications
%        s_amp_thresh = quantile(raw_copy_num_mat(find(diploid_snps_vec==1),i), amp_quantile);
        s_amp_thresh = quantile(raw_copy_num_vec_smooth(find(diploid_snps_vec==1)), amp_quantile);
        if(isnan(s_amp_thresh))
            s_amp_thresh = thresh_amp;
        end
        if (strcmpi(char(gender{i}),'M'))
%             zero_one_mat(chr_x_ind(raw_copy_num_mat(chr_x_ind,i)>s_amp_thresh-two_copy_on_copy_diff),i) = 1;
%             zero_one_mat(non_chr_x_ind(raw_copy_num_mat(non_chr_x_ind,i)>s_amp_thresh),i) = 1;
            zero_one_mat(chr_x_ind(raw_copy_num_vec_smooth(chr_x_ind)>s_amp_thresh-two_copy_on_copy_diff),i) = 1;
            zero_one_mat(non_chr_x_ind(raw_copy_num_vec_smooth(non_chr_x_ind)>s_amp_thresh),i) = 1;
        else
%            zero_one_mat(raw_copy_num_mat(:,i)>s_amp_thresh,i) = 1;
            zero_one_mat(raw_copy_num_vec_smooth>s_amp_thresh,i) = 1;
        end
    end
end

