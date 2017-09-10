%function std_mat = calc_std_mat(copy_num_mat)
function std_mat = calc_std_mat(copy_num_mat)

smooth_param = 30;
% smooth all samples
for i = 1:size(copy_num_mat,2)
    copy_num_mat(:,i) = smooth(copy_num_mat(:,i),smooth_param);
end
num_snps = size(copy_num_mat, 1);
% calc std in moving window
std_frac = 1000;
num_snps_std = floor(num_snps/std_frac);
cont = 1;
ind = 1;
while cont
    start_ind = (ind-1)*num_snps_std+1;
    end_ind = ind*num_snps_std;
    if(end_ind < num_snps)
        std_mat(ind, :) = std(copy_num_mat(start_ind:end_ind, :));
        ind = ind+1;
    else
        cont=0;
    end
end