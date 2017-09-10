qatar_dir = '../../common_disease_model/data/qatar/'; %    qatar_homsnps.txt
qatar_freq_file = fullfile(qatar_dir, 'qatar_homsnps.mat'); % Parse Qatar data: cumulative homozygoucy curves
europe_freq_file = fullfile(qatar_dir, 'caucasian_homsnps.mat');

if(exist(qatar_freq_file, 'file'))
    load(qatar_freq_file);
else
    Q = load(file_name_to_txt(qatar_freq_file));
    save(qatar_freq_file, 'Q');
end
if(exist(europe_freq_file, 'file'))
    load(europe_freq_file);
else
    E = load(file_name_to_txt(europe_freq_file));
    save(europe_freq_file, 'E');
end

% figure; hold on; plot(E); plot(Q, 'r');
% legend('Europe', 'Qatar');
mean_E = mean(E); std_E = std(E); 
mean_Q = mean(Q); std_Q = std(Q); 

legend_vec = []; 
legend_vec{1} = ['Europe, \mu_E = ' num2str(mean_E, 3) ', \sigma_E = ' num2str(std_E, 3)];
legend_vec{2} = ['Qatar, \mu_Q = ' num2str(mean_Q, 3) ', \sigma_Q = ' num2str(std_Q, 3)];

num_bins = length(E);
figure; hold on; plot(sort(E), (1:length(E)) ./ num_bins, 'linewidth', 2);  % plot cumulative 
plot(sort(Q), (1:length(Q)) ./ num_bins, 'r', 'linewidth', 2);
legend(legend_vec, 4); title('Length of heterozygosity runs'); 
xlabel('Frac. SNPs in ROH > 1Mb'); ylabel('Cumulative Distribution'); 
my_saveas(gcf, fullfile(qatar_dir, 'figs', 'qatar_cumulative'), {'epsc', 'fig', 'jpg', 'pdf'});


figure; hold on; % plot density 
hist_density(E, 20, 'b', [], [], 2);
hist_density(Q, 20, 'r', [], [], 2);
hist_density(squareform(IBD_mat), 100, 'g', [], [], 2);
sigma_M = std(squareform(IBD_mat))
mean_M = mean(squareform(IBD_mat))
legend_vec{3} = ['model, \mu_M=' num2str(mean_M) ', \sigma_M=' num2str(sigma_M, 3)];
legend([legend_vec], 1); title('Length of heterozygosity runs'); 
xlabel('Frac. SNPs in ROH > 1Mb'); ylabel('Density distribution'); 
line([mean_Q mean_Q], [0 10], 'linewidth', 2, 'color', 'k');
xlim([0,0.25]);
my_saveas(gcf, fullfile(qatar_dir, 'figs', 'qatar_density'), {'epsc', 'fig', 'jpg', 'pdf'});

figure; hold on; % plot density of only Qatar and model (without europe)
hist_density(Q, 20, 'r', [], [], 2);
hist_density(squareform(IBD_mat), 100, 'g', [], [], 2);
sigma_M = std(squareform(IBD_mat))
mean_M = mean(squareform(IBD_mat))
legend_vec{3} = ['model, \mu_M=' num2str(mean_M) ', \sigma_M=' num2str(sigma_M, 3)];
legend([legend_vec(2:3)], 1); title('Length of heterozygosity runs'); 
xlabel('Frac. SNPs in ROH > 1Mb'); ylabel('Density distribution'); 
line([mean_Q mean_Q], [0 10], 'linewidth', 2, 'color', 'k');
xlim([0,0.25]);
my_saveas(gcf, fullfile(qatar_dir, 'figs', 'qatar_only_density'), {'epsc', 'fig', 'jpg', 'pdf'});

