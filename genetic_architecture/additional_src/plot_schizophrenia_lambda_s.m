% function plot_schizophrenia_lambda_s()
figs_dir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/schizophrenia/';
kinship_vec = [1 0.5 0.25 0.125]; 
familial_vec = {'MZ-twins', 'siblings', 'grand-children', '1st-cousins'}; 

lambda_R_table = [52.1 37.2 44.1; ...
    8.6 8.6 8.6; ...
    3.3 3.36 3.3; ...
    1.8 1.91 1.98];

prevalence = 0.0085; 

for i=1:4
    for j=1:3
        r_R_table(i,j) = kinship_vec(i) * ...
            familial_risk_to_heritability(lambda_R_table(i,j), 'liability', prevalence, kinship_vec(i));
    end
end

legend_vec = {'A_{\Delta}', 'LP_{\Delta}(2)', 'data'}; 

figure; hold on;
plot(kinship_vec, log2(lambda_R_table(:,2)), 'r', 'linewidth', 2); % additive (liability-threshold) model
plot(kinship_vec, log2(lambda_R_table(:,3)), 'g', 'linewidth', 2); % LP(2) model
plot(kinship_vec, log2(lambda_R_table(:,1)), '*', 'markersize', 8); % data
legend(legend_vec, 2); legend('boxoff');
xlabel('$2 \times \mathrm{kinship} \: (\kappa_R)$', 'interpreter', 'latex', 'fontsize', 15); 
ylabel('$\mathrm{Familial \: \; risk} \: \log(\lambda_R)$', 'interpreter', 'latex', 'fontsize', 15);
title([repmat(' ', 1, 80) 'supp-fig7a'], 'fontsize', 16, 'fontweight', 'bold');
xlim([0,1]);
x_vec = kinship_vec+0.03; x_vec(1) = x_vec(1)-0.2;
y_vec = log2(lambda_R_table(:,1))-0.01;
for i=1:4 % add text
    text(x_vec(i), y_vec(i), familial_vec{i}, 'fontsize', 14); 
end

my_saveas(gcf, fullfile(figs_dir, 'schizophrenia_familial_risk_fit'), {'epsc', 'jpg', 'fig'});



figure; hold on; % plot in correlation (liability) space 
plot(kinship_vec, r_R_table(:,2), 'r', 'linewidth', 2); % additive (liability-threshold) model
plot(kinship_vec, r_R_table(:,3), 'g', 'linewidth', 2); % LP(2) model
plot(kinship_vec, r_R_table(:,1), '*', 'markersize', 8); % data
% legend(legend_vec, 2);  legend('boxoff'); % no need for legend again
xlabel('$2 \times \mathrm{kinship} \: (\kappa_R)$', 'interpreter', 'latex', 'fontsize', 15); 
ylabel('$\mathrm{Phenotypic \:\: correlation} \: (r_R)$', 'interpreter', 'latex', 'fontsize', 15);
title([repmat(' ', 1, 80) 'supp-fig7b'], 'fontsize', 16, 'fontweight', 'bold');
xlim([0,1]); ylim([0,1]);
y_vec = r_R_table(:,1)-0.01;
for i=1:4 % add text
    text(x_vec(i), y_vec(i), familial_vec{i}, 'fontsize', 14); 
end
my_saveas(gcf, fullfile(figs_dir, 'schizophrenia_familial_correlation_fit'), {'epsc', 'jpg', 'fig'});




