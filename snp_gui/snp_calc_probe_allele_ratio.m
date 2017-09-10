%function [ratio_vec] = snp_calc_probe_allele_ratio(probe_intens_mat, proj_name, mean_mm_flag)
function [ratio_vec] = snp_calc_probe_allele_ratio(probe_intens_mat, mean_mm_flag) 

[pa_vec, ma_vec, pb_vec, mb_vec] = snp_get_probe_pm_mm_columns(proj_name);
                        get_probe_pm_mm_ind(column_content_cell, quartet_num_vec, match_cell, strand_cell);

%Calculate PM-MM for each probe pair or probe - MM average:
%d = (PM(i)-MM(i)) or (PM(i)-avg(MM))

num_snp = size(probe_intens_mat, 1);
if(mean_mm_flag)
    mm_avg_mat = (probe_intens_mat(:, ma_vec)+probe_intens_mat(:, mb_vec))/2;
    d_mat_a = probe_intens_mat(:, pa_vec)- mm_avg_mat;
    d_mat_b = probe_intens_mat(:, pb_vec)- mm_avg_mat;
else
    d_mat_a = (probe_intens_mat(:, pa_vec)-probe_intens_mat(:, ma_vec));
    d_mat_b = (probe_intens_mat(:, pb_vec)-probe_intens_mat(:, mb_vec));
end

probe_intens_mat_filtered = probe_intens_mat;
probe_intens_mat_filtered_pa = probe_intens_mat_filtered(:, pa_vec);
probe_intens_mat_filtered_pb = probe_intens_mat_filtered(:, pb_vec);
probe_intens_mat_filtered_pa(find(d_mat_a<0)) = nan;
probe_intens_mat_filtered_pb(find(d_mat_b<0)) = nan;
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
% for i = 1:num_quartets
%     %quartets_avg_mm(:,i) = mean(probe_intens_mat_filtered(:, [ma_vec(i) mb_vec(i)])')';
%     if(mean_mm_flag)
%         quartets_avg_mm(:,i) = mean_not_nan(probe_intens_mat_filtered(:, [ma_vec(i) mb_vec(i)])')';
%         quartets_a(:,i) = probe_intens_mat_filtered(:, pa_vec(i))-quartets_avg_mm(:,i);
%         quartets_b(:,i) = probe_intens_mat_filtered(:, pb_vec(i))-quartets_avg_mm(:,i);
%     else
%         quartets_a(:,i) = probe_intens_mat_filtered(:, pa_vec(i))-probe_intens_mat_filtered(:, ma_vec(i));
%         quartets_b(:,i) = probe_intens_mat_filtered(:, pb_vec(i))-probe_intens_mat_filtered(:, mb_vec(i));
%     end
% end

a_intens = mean_not_nan(quartets_a')';
b_intens = mean_not_nan(quartets_b')';
%ratio_vec = a_intens./(a_intens+b_intens);
ratio_vec = a_intens./b_intens;