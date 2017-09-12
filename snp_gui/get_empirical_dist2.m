function ret = get_empirical_dist2(empirical_dist_data, observed_data_2)

ret = 0;

num_emp_obs = length(empirical_dist_data);
empirical_dist_data_pairs = zeros(num_emp_obs-1,2);
if(size(empirical_dist_data,1) ==1)
    empirical_dist_data = empirical_dist_data';
end
empirical_dist_data_pairs(:,1) = empirical_dist_data(1:end-1);
empirical_dist_data_pairs(:,2) = empirical_dist_data(2:end);

ret = length(intersect(find(empirical_dist_data_pairs(:,1) <= observed_data_2(1)), ...
    find(empirical_dist_data_pairs(:,2) <= observed_data_2(2))));

ret = ret / (num_emp_obs-2+1);
