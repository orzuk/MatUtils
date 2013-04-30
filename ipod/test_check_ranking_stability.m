% Generate graphs according to different models.
% and compute the stability of them to sampling effects

SAMP_EDGES = 0; SAMP_NODES = 1;  % sampling flags

ttt=cputime;
alpha_vec = 0.01:0.01:0.2;  % enumerate different values of alpha


N=100; alpha=0.3; iters =500; p_vec1 = [0:0.01:0.25];  samp_type = SAMP_EDGES; m=4; p_vec2 = (2*m / (N-2*m)) .* p_vec1; % make almost symmetric noise


f_mean_mat_ER = zeros(length(alpha_vec), length(p_vec));
f_mean_mat_BA = zeros(length(alpha_vec), length(p_vec));

i=1;
for alpha = 0.3% alpha_vec
    g_t=0; [f_mat_ER corr_mat_ER] = check_ranking_stability(N, alpha, iters, p_vec, p_vec2, samp_type, g_t, 2*m/N);  % Erdos-Renyi
    g_t=2; [f_mat_BA corr_mat_BA] = check_ranking_stability(N, alpha, iters, p_vec, p_vec2, samp_type, g_t, m);  % Barabasi-Albert

    
%     h = hist(r, max(r)-min(r)+1);
%     figure;  xlabel('k'); ylabel('P(k)'); hold on; loglog(min(r):max(r), h./N, '*');   
    
    f_mean_mat_ER(i,:) = mean(f_mat_ER,2);
    f_mean_mat_BA(i,:) = mean(f_mat_BA,2);

    i=i+1;
end



% Do plot
figure; hold on; title(['overlap in top degrees \alpha=' num2str(alpha) ' N = ' num2str(N)]); xlabel('p'); ylabel('f');
plot(p_vec, mean(f_mat_ER,2)); plot(p_vec, mean(f_mat_BA,2), 'r');
plot(p_vec, mean(corr_mat_ER,2), 'g'); plot(p_vec, mean(corr_mat_BA,2), 'c');
legend('ER', 'BA', 'ER corr', 'BA corr');


ttt = cputime - ttt

return; % Stop

ind=500; p=p_vec(ind);
figure; hold on; title(['overlap in top degrees p=' num2str(p) ' N = ' num2str(N)]); xlabel('\alpha'); ylabel('f');
plot(alpha_vec, f_mean_mat_ER(:,ind)); plot(alpha_vec, f_mean_mat_BA(:,ind), 'r'); legend('ER', 'BA');


figure; subplot(2,1,1);  imagesc(f_mean_mat_ER); colorbar; xlabel('p'); ylabel('\alpha'); title('Erdos-Renyi');
subplot(2,1,2);  imagesc(f_mean_mat_BA); colorbar; xlabel('p'); ylabel('\alpha'); title('Barabasi-Albert');


