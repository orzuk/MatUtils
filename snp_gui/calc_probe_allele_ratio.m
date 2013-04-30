%function [ratio_vec] = calc_probe_allele_ratio(probe_intens_mat, proj_name, mean_mm_flag)
function [ratio_vec] = calc_probe_allele_ratio(probe_intens_mat, mean_mm_flag, ...
                       pa_vec, ma_vec, pb_vec, mb_vec) 

EPSILON = 0.0000000000001; % Used to avoid infinite and zero ratios                   
                   
%Calculate PM-MM for each probe pair or probe - MM average:
%d = (PM(i)-MM(i)) or (PM(i)-avg(MM))

num_snp = size(probe_intens_mat, 1);
if(mean_mm_flag)
%    mm_avg_mat = (probe_intens_mat(:, ma_vec)+probe_intens_mat(:, mb_vec))/2;
    mm_avg_mat = add_non_nan_ind(probe_intens_mat(:, ma_vec),probe_intens_mat(:, mb_vec)); % this function adds only the non nan indices AND calculate average
    d_mat_a = probe_intens_mat(:, pa_vec)- mm_avg_mat;
    d_mat_b = probe_intens_mat(:, pb_vec)- mm_avg_mat;
else
    d_mat_a = (probe_intens_mat(:, pa_vec)-probe_intens_mat(:, ma_vec));
    d_mat_b = (probe_intens_mat(:, pb_vec)-probe_intens_mat(:, mb_vec));
end

probe_intens_mat_filtered = probe_intens_mat;
probe_intens_mat_filtered_pa = probe_intens_mat_filtered(:, pa_vec);
probe_intens_mat_filtered_pb = probe_intens_mat_filtered(:, pb_vec);
probe_intens_mat_filtered_pa(find(d_mat_a<=0)) = nan;
probe_intens_mat_filtered_pb(find(d_mat_b<=0)) = nan;
probe_intens_mat_filtered(:, pa_vec) = probe_intens_mat_filtered_pa;
probe_intens_mat_filtered(:, pb_vec) = probe_intens_mat_filtered_pb;

%The mismatch average for each quartet is calculated, MMavg.
if(mean_mm_flag)
    quartets_a = probe_intens_mat_filtered(:, pa_vec) - mm_avg_mat;
    quartets_b = probe_intens_mat_filtered(:, pb_vec) - mm_avg_mat;
else
    quartets_a = probe_intens_mat_filtered(:, pa_vec) - probe_intens_mat_filtered(:, ma_vec);
    quartets_b = probe_intens_mat_filtered(:, pb_vec) - probe_intens_mat_filtered(:, mb_vec);
end
 
a_intens = nanmean(quartets_a')';
b_intens = nanmean(quartets_b')';
a_intens(isnan(a_intens)) = 0;
b_intens(isnan(b_intens)) = 0;
%ratio_vec = a_intens./(a_intens+b_intens);
ratio_vec = single( (a_intens+EPSILON)./ (b_intens+EPSILON) ); % transfer to singles to save space
