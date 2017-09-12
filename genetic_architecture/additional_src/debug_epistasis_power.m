% function debug_epistasis_power()
h_x=0.1; % This is the heritability of each of two loci for the LT model for the liability.
mu=0.5; N=2;
x_mu = norminv(1-mu);
[M S S_ONE] = z_stats_given_two_x_MLT(N,1,mu, h_x), 
[M_LT S_LT] = z_stats_given_two_x_MLT(1,1,mu, h_x), 
TOTAL_EFFECT = 100*(M*(1-M)-S)/(M*(1-M))
MAIN_EFFECT = 100*(M*(1-M)-S_ONE)*2/(M*(1-M))
INTERACTION = 100*(S_LT-S)/(M*(1-M))
MAIN_EFFECT_LT = 100*(M*(1-M)-S_LT)/(M*(1-M))
n_samples = 10000;

XZ = simulate_genotype_phenotype( [0.5 sqrt(h_x) 1 mu], ...
    n_samples, n_samples, ...
    1, 'binary', 'epistasis', 1);% Simulate data
h_x_MLT = N*heritability_scale_change_MLT(h_x,1,N,mu,'MLT'); % The heritability of each of two liabilities. we use two x1,x2 on the same liability
XZ2 = simulate_genotype_phenotype( [0.5 sqrt(h_x_MLT) N mu], ...
    n_samples, n_samples, ...
    1, 'binary', 'epistasis', 1);% Simulate data

C = corr(XZ(1,:,1)', XZ(1,:,3)'), C_again = corr(XZ(1,:,2)', XZ(1,:,3)')
C2 = corr(XZ2(1,:,1)', XZ2(1,:,3)'), C2_again = corr(XZ2(1,:,2)', XZ2(1,:,3)')

figure; hold on; bins2d = linspace(min(XZ(1,:,1))*1.1, max(XZ(1,:,1))*1.1, 25); 
disease_inds = find(XZ(1,:,3));
total_hist = hist2d(reshape(XZ(1,:,[1 2]), length(XZ), 2), bins2d, bins2d) 
diseased_hist = hist2d(reshape(XZ(1,disease_inds,[1 2]), length(disease_inds), 2), bins2d, bins2d); 
Plot2dHist(diseased_hist ./ max(1, total_hist), bins2d, bins2d, 'X_1', 'Y_1', ...
    ['Disease state under LT. \mu=' num2str(mu*100,3) '%, h_x^2=' num2str(h_x*100,3) '%']);

figure; hold on; bins2d = linspace(min(XZ2(1,:,1))*1.1, max(XZ2(1,:,1))*1.1, 25); 
disease_inds = find(XZ2(1,:,3));
total_hist = hist2d(reshape(XZ2(1,:,[1 2]), length(XZ2), 2), bins2d, bins2d) 
diseased_hist = hist2d(reshape(XZ2(1,disease_inds,[1 2]), length(disease_inds), 2), bins2d, bins2d); 
Plot2dHist(diseased_hist ./ max(1, total_hist), bins2d, bins2d, 'X_1', 'Y_1', ...
    ['Disease state under MLT(2,1) \mu=' num2str(mu*100,3) '%, h_x^2=' num2str(h_x*100,3) '%']);


figure; hold on; 
scatter(XZ(1,:,1), XZ(1,:,2), 2.5, XZ(1,:,4)); colorbar;
title(['Disease prob. under LT. \mu=' num2str(mu*100,3) '%, h_x^2=' num2str(h_x*100,3) '%']);
xlabel('X_1'); ylabel('X_2');
xlim([-4 4]); ylim([-4 4]);

figure; hold on; 
scatter(XZ2(1,:,1), XZ2(1,:,2), 2.5, XZ2(1,:,4)); colorbar;
title(['Disease prob. under MLT(2,1) \mu=' num2str(mu*100,3) '%, h_x^2=' num2str(h_x*100,3) '%']);
xlabel('X_1'); ylabel('X_2');
xlim([-4 4]); ylim([-4 4]);

figure; hold on; 
XZ2_expected_LT = 1 - (normcdf( (x_mu - XZ2(:,:,1) .* sqrt(h_x) - ...
                            XZ2(:,:,2) .* sqrt(h_x)) ./ ...
                        sqrt(1-2*h_x))); % z_expected
scatter(XZ2(1,:,1), XZ2(1,:,2), 2.5, XZ2(1,:,4)-XZ2_expected_LT); colorbar;
title(['Disease prob. diff MLT(2,1)-LT \mu=' num2str(mu*100,3) '%, h_x^2=' num2str(h_x*100,3) '%']);
xlabel('X_1'); ylabel('X_2');
xlim([-4 4]); ylim([-4 4]);


% plot(XZ(1,disease_inds,1), XZ(1,disease_inds,2), 'r.'); 
% title(['Disease state under LT. \mu=' num2str(mu*100,3) '%, h_x^2=' num2str(h_x*100,3) '%']);
% xlabel('X_1'); ylabel('X_2'); % legend('healthy', 'diseased'); 
% 
