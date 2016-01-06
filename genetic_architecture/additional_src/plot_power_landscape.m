% Plot detection power for different effect sizes
%   Detailed explanation goes here

function [power_mat fig_handle] = plot_power_landscape(num_cases, num_controls, ...
    alpha, trait_type, fig_outfile)

AssignGeneralConstants;
x_vec = 0.01:0.002:0.5; % MAF
beta_vec = 0:0.0005:1; % effect size

power_mat = zeros(length(x_vec), length(beta_vec));

test_type = 'single-locus'; test_stat = 'chi-square-QTL-analytic';
for i=1:length(x_vec) % loop on MAF 
    run_power_i = i
    p_z_x_mat = QTL_params_to_p_mat(repmat(x_vec(i), length(beta_vec), 1), beta_vec', 1);
    power_mat(i,:) = compute_association_power(p_z_x_mat, num_cases, num_controls, ...
        alpha, [], test_type, test_stat, 'population');
end

figure; imagesc_with_labels(power_mat', x_vec, beta_vec, [], [], [25 250]); colormap('gray');
xlabel('MAF'); ylabel('\beta'); 
title(['QTL Detection Power for n = ' num2str(num_cases) ' with \alpha = ' ...
    num2str(alpha) ' significance threshold']);  
my_saveas(gcf, fig_outfile, format_fig_vec); 
fig_handle = []; 
