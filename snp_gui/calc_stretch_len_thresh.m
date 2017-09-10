% Compute stretch width threshold 
%function stretch_len_thresh = calc_stretch_len_thresh(zero_one_vec, joined_ones_vec, frac_cut, max_stretch_len)
function stretch_len_thresh = calc_stretch_len_thresh(zero_one_vec, joined_ones_vec, frac_cut, max_stretch_len)

start_num_joined = 1;
max_stretch_len = 1000;
num_aberr = sum(zero_one_vec);

rand_consequent_vec = get_rand_consequent_vec(length(zero_one_vec), num_aberr, start_num_joined, max_stretch_len);
joined_ones_vec(find(joined_ones_vec < start_num_joined)) = [];
stretches_hist_rand = histc(rand_consequent_vec, [start_num_joined: max_stretch_len]);
stretches_hist = histc(joined_ones_vec, [start_num_joined: max_stretch_len]);
% find stretch len such that the rest of rand fraction <= frac_cut
rand_frac_vec = zeros(1, length(stretches_hist_rand));

rand_rel_real_frac_vec = zeros(1, length(stretches_hist_rand));
num_rand_stretches = sum(stretches_hist_rand);
max_rand_stretch = max(find(stretches_hist_rand>0));
for i = 1:max_rand_stretch
    rand_frac_vec(i) = sum(stretches_hist_rand(i:end))/num_rand_stretches;
    if(sum(stretches_hist(i:end)) ~=0)
        rand_rel_real_frac_vec(i) = sum(stretches_hist_rand(i:end))/sum(stretches_hist(i:end));
    end
end

%ind = find(rand_frac_vec <= frac_cut);
ind = find(rand_rel_real_frac_vec <= frac_cut);
if(length(ind)>0)
   stretch_len_thresh = start_num_joined+min(ind)-1;
else
    stretch_len_thresh = max_stretch_len;
end