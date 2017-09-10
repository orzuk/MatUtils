% Here we load a data, compute the correlations and get the desired
% fraction
path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
path(path, 'E:\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric\nmm\linalg');


TOL = 0.00000000001;

% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;
rand_flag = GAUSSIAN;

% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TRUE_AND_SAMPLED;

max_nsamples = 200; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)

res = 20; % the resulution of N that we use
alpha = 0.014; alpha_vec =  [0.012  0.12]; % corresponds to ~70 and ~700 genes

%nsamples=3;
iters = 1200;

nsamples_vec = res:res:max_nsamples;


load('../data/breastData.mat'); % The name of the data file containing the expression and labels.


% Now call the function to do all work
[N_VEC Pearson_corrs Fisher_Zs Rho_sig_mean Rho_sig_std Rho_bias_mean Rho_bias_std Z_sig_mean ...
    Z_sig_std Z_bias_mean Z_bias_std chunk_Rho_sig chunk_Z_sig chunk_Rho_bias chunk_Z_bias hist_Pearson_Ps hist_Fisher_Zs random_genes_picked] = ...
    CalcZVarianceFromDataFunc(R, alpha, rand_flag, iters);


% Now plot results
figure; subplot(2,2,1);  hold on; title('Pearson P std, from Breast data ');
errorbar(N_VEC, Rho_sig_mean(N_VEC),Rho_sig_std(N_VEC));
xlabel('N'); ylabel('P std'); %legend('from data', 'from Fisher');


subplot(2,2,2);  hold on; title('Fisher Z std, from Breast data and analytic');
plot(N_VEC, 1./sqrt(N_VEC-3), 'r');
errorbar(N_VEC, Z_sig_mean(N_VEC),Z_sig_std(N_VEC));

xlabel('N'); ylabel('Z std'); legend( '(N-3)^{-1/2}','from data');

subplot(2,2,3);  hold on; title('Pearson P bias, from Breast data');
errorbar(N_VEC, Rho_bias_mean(N_VEC),Rho_bias_std(N_VEC));
xlabel('N'); ylabel('P bias'); %legend('from data', 'from Fisher');


subplot(2,2,4);  hold on; title('Fisher Z bias, from Breast data');
errorbar(N_VEC, Z_bias_mean(N_VEC),Z_bias_std(N_VEC));
xlabel('N'); ylabel('Z bias'); %%legend('from data', 'from Fisher');



% Plot the same thing without the error bars
figure; subplot(2,2,1);  hold on; title('Pearson P std, from Breast data ');
plot(N_VEC, Rho_sig_mean(N_VEC));
xlabel('N'); ylabel('P std'); %legend('from data', 'from Fisher');


subplot(2,2,2);  hold on; title('Fisher Z std, from Breast data and analytic');
plot(N_VEC, 1./sqrt(N_VEC-3), 'r');
plot(N_VEC, Z_sig_mean(N_VEC));

xlabel('N'); ylabel('Z std'); legend( '(N-3)^{-1/2}','from data');

subplot(2,2,3);  hold on; title('Pearson P bias, from Breast data');
plot(N_VEC, Rho_bias_mean(N_VEC));
xlabel('N'); ylabel('P bias'); %legend('from data', 'from Fisher');


subplot(2,2,4);  hold on; title('Fisher Z bias, from Breast data');
plot(N_VEC, Z_bias_mean(N_VEC));
xlabel('N'); ylabel('Z bias'); %%legend('from data', 'from Fisher');




% Plot to see correlation between the 'true' value and the variance
figure; subplot(2,2,1); hold on; plot(Pearson_corrs, chunk_Rho_sig, '.'); title(['True Pearson P and estimator std. N=' num2str(N_VEC(end))]);
xlabel('true P'); ylabel('std');
subplot(2,2,2); hold on; plot(Fisher_Zs, chunk_Z_sig, '.'); title(['True Fisher Z and estimator std. N='  num2str(N_VEC(end))]);
xlabel('true Z'); ylabel('std');


subplot(2,2,3); hold on; plot(Pearson_corrs, chunk_Rho_bias, '.'); title(['True Pearson P and estimator bias.  N=' num2str(N_VEC(end))]);
xlabel('true P'); ylabel('bias');
subplot(2,2,4); hold on; plot(Fisher_Zs, chunk_Z_bias, '.'); title(['True Fisher Z and estimator bias. N='  num2str(N_VEC(end))]);
xlabel('true Z'); ylabel('bias');


% Now plot the histograms !!!
figure; subplot(4,4,1);

for i=1:4
    for j=1:4
        subplot(4,4,4*(i-1)+j); hold on;
        hist(hist_Pearson_Ps(4*(i-1)+j,:), 100); title(['Pearson hist gene # ' num2str(random_genes_picked(4*(i-1)+j))]);
    end
end


figure; subplot(4,4,1);

for i=1:4
    for j=1:4
        subplot(4,4,4*(i-1)+j); hold on;
        hist(hist_Fisher_Zs(4*(i-1)+j,:), 100); title(['Fisher hist gene # ' num2str(random_genes_picked(4*(i-1)+j))]);
    end
end





