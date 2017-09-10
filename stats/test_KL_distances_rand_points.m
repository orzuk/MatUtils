% Script for runnning the test_KL_distances_rand_points function


ttt = cputime;
plot_results = 0;
max_N = 30; iters = 2;
% % % % ttt = cputime;
% % % % KL_mat_rand = KL_distances_rand_points(max_N, iters,0) .* log(2); % Two random points. Here we can do more
% % % % rand_t = cputime - ttt
% % % % ttt = cputime;
% % % % KL_mat_half_rand = KL_distances_rand_points(max_N, iters,0.5)  .* log(2); %  A random point and the closest low-dim point
% % % % half_rand_t = cputime - ttt
% % % % % ttt = cputime;
% % % % % KL_mat_Cij = KL_distances_rand_points(max_N, iters,1) .* log(2);
% % % % % Cij_t = cputime - ttt
% % % % ttt = cputime;
% % % % KL_mat_max_k = KL_distances_rand_points(max_N, iters,2); % .* log(2); %  A low-dim random point anda random low-dim point
% % % % KL_mat_max_k_closest = KL_mat_max_k{1} .* log(2); KL_mat_max_k = KL_mat_max_k{2} .* log(2);
% % % % rand_max_k_t = cputime - ttt


% ttt = cputime;
% KL_mat_max_k_closest = KL_distances_rand_points(max_N, iters,2.5) .* log(2); %  A low-dim random point and a random low-dim point
% closest_max_k_t = cputime - ttt
ttt = cputime;
KL_mat_max_k_low_model = KL_distances_rand_points(max_N, iters,3); % A random point and a random low-dim point
KL_mat_max_k_low_model_closest = KL_mat_max_k_low_model{2} .* log(2);
KL_mat_max_k_low_model_no_rand = KL_mat_max_k_low_model{3} .* log(2);
KL_mat_max_k_low_model = KL_mat_max_k_low_model{1} .* log(2);
max_k_low_model_t = cputime - ttt


% Do not save, not to 'ruin' the bigg run
%save('AllKLMatsTrees.mat', 'KL_mat_rand', 'KL_mat_half_rand', 'KL_mat_max_k', 'KL_mat_max_k_closest', 'max_k_low_model_t', 'KL_mat_max_k_low_model_closest', 'KL_mat_max_k_low_model_no_rand');


if(plot_results) % plot results:
    figure; hold on;
    errorbar(1:max_N, mean(KL_mat_rand'), std(KL_mat_rand'));
    errorbar([1:max_N]+0.05, mean(KL_mat_half_rand'), std(KL_mat_half_rand'),'--');
    errorbar([1:max_N]+0.1, mean(KL_mat_max_k'), std(KL_mat_max_k'),'g');
    errorbar([1:max_N]+0.15, mean(KL_mat_max_k_closest'), std(KL_mat_max_k_closest'),'c');
    errorbar([1:max_N]+0.2, mean(KL_mat_max_k_low_model'), std(KL_mat_max_k_low_model'),'m');
    errorbar([1:max_N]+0.2, mean(KL_mat_max_k_low_model_closest'), std(KL_mat_max_k_low_model_closest'),'r');
    %%errorbar([1:max_N]+0.3, mean(KL_mat_Cij'), std(KL_mat_Cij'),'r');
    legend('Two rands', 'Half Rand', 'Random InDegree 1',  'Closest InDegree 1',  'Random Indegree 1 low model', 'Closest Indegree 1 low model');  %%5legend('random', 'Cij', 'InDegree 2');
    title('KL distance from a random point'); xlabel('n'); ylabel('Relative Entropy (NATS)');
    
    % Show harmonic series for comparison
    % har_vec = [0 cumsum(1./[2:2^max_N])]; har_vec  = har_vec(2.^[1:max_N])./log(2);
    % plot(1:max_N, har_vec, 'k');
    %
    % elapsed_time = cputime-ttt
    
    plot(1:max_N, ones(1,max_N), 'k:'); plot(1:max_N, ones(1,max_N)-0.577215665, 'k:');
    
    % Plot only the low-low dim errors
    figure; hold on; plot([3:max_N], log2(mean(KL_mat_max_k_low_model(3:max_N,:)')),'m');
    plot([3:max_N], log2(mean(KL_mat_max_k_low_model_closest(3:max_N,:)')),'r');
    plot([3:max_N], log2(mean(KL_mat_max_k_low_model_no_rand(3:max_N,:)')),'b');
    title('KL distnace from random low dimensional point'); xlabel('n'); ylabel('log(Relative Entropy)');
    legend('Random Indegree 1 low model', 'Closest Indegree 1 low model', 'Random without Rand');
    
end

