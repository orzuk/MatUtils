%function [whole_chr_aberr_flag del_flag amp_flag] = check_if_whole_chr_aberr(chr_copy_num_vec, del_amp_flag)
function [whole_chr_aberr_flag del_flag amp_flag] = check_if_whole_chr_aberr(chr_copy_num_vec, del_amp_flag, chr_num, gender, frac_whole_arm)

%frac_whole_arm = 0.3;
whole_chr_aberr_flag = 0;
del_flag = 0;
amp_flag = 0;

if strcmpi(gender,'M') && chr_num==23
    thresh=1;
else
    thresh=2;
end

if(length(find(chr_copy_num_vec < thresh)) > frac_whole_arm*length(chr_copy_num_vec))
    del_flag = 1;
end

if(length(find(chr_copy_num_vec > thresh)) > frac_whole_arm*length(chr_copy_num_vec))
    amp_flag = 1;
end

if(del_amp_flag==1 & del_flag==1) % deletion
    whole_chr_aberr_flag = 1;
end

if(del_amp_flag==2 & amp_flag==1) % amplification
    whole_chr_aberr_flag = 1;
end
