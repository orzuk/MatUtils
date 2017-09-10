%function [ret, aberr_str] = check_if_chr_loc_big_aberr(big_aberr_mat, end_p_location, chr_num, chr_loc_vec);
function [ret, aberr_str] = check_if_chr_loc_big_aberr(big_aberr_mat, end_p_location, chr_num, chr_loc_vec);

aberr_str = {'no aberr'};

ret = 0;
% check if there is a big aberration in this chromosome
if(sum(abs(big_aberr_mat(chr_num,:)))>0)
    % check if chr_loc_vec in p\q\both
    p_flag = length(find(chr_loc_vec<= end_p_location(chr_num))>0);
    q_flag = length(find(chr_loc_vec> end_p_location(chr_num))>0);
    if(p_flag)
        if(big_aberr_mat(chr_num,1)==1) % amplification
            aberr_str = 'AMP';
            ret = 1;
        end
        if(big_aberr_mat(chr_num,1)==-1) % deletion
            aberr_str = 'DEL';
            ret = 1;
        end
    end
    if(q_flag)
        if(big_aberr_mat(chr_num,2)==1) % amplification
            aberr_str = 'AMP';
            ret = 1;
        end
        if(big_aberr_mat(chr_num,2)==-1) % deletion
            aberr_str = 'DEL';
            ret = 1;
        end
    end
end
