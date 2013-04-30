% The function calculate some statistics for random instances
% of the graph alignment problem.
function [local_min_num_mean local_min_num_std global_min_num_mean global_min_num_std] = statistics_for_ga(N_vec,num_graphs, graph_type, param)


global_min_num = zeros(length(N_vec), num_graphs);
local_min_num = zeros(length(N_vec), num_graphs);
ground_state_energies = zeros(length(N_vec), num_graphs);

N_i=1;
for N=N_vec
    NN=N
    for g=1:num_graphs
        G1 = generate_random_graph(graph_type, N, param);
        G2 = generate_random_graph(graph_type, N, param);
        [global_min_num(N_i, g) local_min_num(N_i, g) ground_state_energies(N_i,g)] =  num_ga_local_minima(G1,G2);
    end
    N_i=N_i+1;
end



global_min_num_mean = mean(global_min_num,2); global_min_num_std = std(global_min_num,[],2);
local_min_num_mean = mean(local_min_num,2); local_min_num_std = std(local_min_num,[],2);
ground_state_energies_mean = mean(ground_state_energies,2); ground_state_energies_std = std(ground_state_energies,[],2);

figure; subplot(3,1,1); hold on;errorbar(N_vec,  global_min_num_mean, global_min_num_std);title('Global Minima'); xlabel('N');
subplot(3,1,2); hold on;errorbar(N_vec,  local_min_num_mean, local_min_num_std);title('Local Minima'); xlabel('N');
subplot(3,1,3); hold on;errorbar(N_vec,  ground_state_energies_mean, ground_state_energies_std);title('Ground State Energy'); xlabel('N');

figure; subplot(length(N_vec),1,1);
for i=1:length(N_vec)
    subplot(length(N_vec),1,i); hold on; hist(ground_state_energies(i,:)); title(['N=' num2str(N_vec(i))]); xlabel('Val.'); ylabel('Prob.');
end



% figure; subplot(2,1,1); hold on;plot(N_vec,  global_min_num_mean);title('Global Minima'); xlabel('N');
% subplot(2,1,2); hold on;plot(N_vec,  local_min_num_mean);title('Local Minima'); xlabel('N');
