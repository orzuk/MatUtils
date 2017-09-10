% Make figures from ESP:
figs_dir = 'C:\Users\user\Dropbox\HUJI\Grants\ISF\2015_16\figs\'

n_genes = 18617;
s_vec = -log2(rand(n_genes, 1));

s_vec = s_vec ./ (5000*max(s_vec));

s_vec(1:1000) = s_vec(1:1000) * 100;

s_vec1 = s_vec + randn(n_genes, 1) .* sqrt(s_vec);
s_vec2 = s_vec + randn(n_genes, 1) .* 0.5 .* sqrt(s_vec);

s_vec1(2000:3000) = s_vec1(2000:3000) * 0.1;
s_vec2(5000:6000) = s_vec2(5000:6000) * 0.2;

s_vec1 = min(max(s_vec1, 0.0000001), 1);
s_vec2 = min(max(s_vec2, 0.0000001), 1);



figure; loglog(s_vec1, s_vec2, '.');
xlabel('Europe', 'fontsize', 16); ylabel('Africa', 'fontsize', 16);
title('Estimated selection coefficients for European and African samples from ESP', 'fontsize', 16);
my_saveas(gcf, fullfile(figs_dir, 'ESP_EuropeAfrica'), {'epsc', 'pdf', 'fig', 'jpg'})






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_plot_file_name = '_s_confidence_interval';
mu = 10^(-5); % rate used for simulations per gene

s_vec = logspace(-4, -1, 1000); 
est_s_vec = s_vec 

z_vec = ones(size(s_vec))
z_vec2 = max(0.0001, s_vec - 0.002-5*s_vec.^2 + 4.*s_vec.^3 + 0.0002 ./ (0.11-s_vec)); 
z_vec2(800:end) = z_vec2(800:end) +  (z_vec2(800:end) - z_vec2(800))*0.5;

ratio_vec = zeros(size(s_vec));
ratio_vec(100:end) = linspace(0, 1, 901)
z_vec_minus = max(0.0001, z_vec2 .* ratio_vec);
ratio_vec2 = ones(size(s_vec));
ratio_vec2(1:end) = linspace(10, 1, 1000)

z_vec_plus = max(0.0000, max(0.0001, z_vec2) .* ratio_vec2);
z_vec_plus(1:149) = 0.75*10^(-6) ./ z_vec_plus(1:149); 

figure; 
loglog(s_vec, z_vec2, 'linewidth', 2); hold on;
loglog(s_vec, s_vec, '--', 'linewidth', 2); hold on;
loglog(s_vec, z_vec_minus, 'k--', 'linewidth', 2);
loglog(s_vec, z_vec_plus, 'k--', 'linewidth', 2);

%                set(gca, 'xscale', 'log'); set(gca, 'yscale', 'log');
xlabel('Estimated selection coefficient $\hat{s}$', 'fontsize', 14, 'fontweight', 'bold', 'interpreter', 'latex');
ylabel('True selection coefficient $s$', 'fontsize', 14, 'fontweight', 'bold', 'interpreter', 'latex');
title(['Distribution of Estimated s'], 'fontsize', 14, 'fontweight', 'bold', 'interpreter', 'latex'); %   Mean=' num2str(mean(f_null)) '.St.d.=' num2str(std(f_null))]);
h_leg = legend({'Median $\hat{s}$', 'True $s$', '$95\%$ confidence'}, 4, 'interpreter', 'latex');
set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
add_faint_grid(0.5);
my_saveas(gcf, fullfile(figs_dir, 'new_s_estimation'), {'epsc', 'pdf', 'fig', 'jpg'})




% Power boost: 
power_vec = 1-2.*log2(rand(n_genes, 1));
power_vec(1:1000) = power_vec(1:1000)*20;
power_vec(5000:10000) = power_vec(5000:10000)+randn(5001,1);

power_vec = max(0.999, power_vec);
figure; hist((power_vec), 1000); 
set(gca,'xscale','log'); hold on; 
xlabel('Power-ratio', 'fontsize', 16); ylabel('Num. genes', 'fontsize', 16); 
title('Relative reduction in sample size due to gene-specific threshold for different genes'); 
my_saveas(gcf, fullfile(figs_dir, 'power_boost'), {'epsc', 'pdf', 'fig', 'jpg'})



