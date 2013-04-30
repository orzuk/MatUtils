% We let X count the number of permutation satysfying 
% that their alignment score is above a certain threshold.
% The we calculate the first two moments of X 
% N is the number of nodes
% p is the edge probability (We take Erdosh-Renyi graphs)
% alpha is the required 'overlap', i.e. how many 'conatrains' we 
% need to satisfy, or, how many edges are correctly aligned
% Enumerate all graphs
function [x_mean x_std p_x_positive_exact p_x_positive_2nd_moment perms_corrs] = calc_x_aligned_moments(N, p, alpha_vec);


num_edges = N*(N-1)/2;


% Now generate all graphs
all_graphs = zeros(2^num_edges, num_edges);
for i=1:num_edges
    all_graphs(:,i) = mod( floor([0:2^num_edges-1] ./ (2^(i-1))), 2);
end


ttt=cputime;

% Calculate the probability of each graph
graphs_probs = zeros(1, 2^num_edges); graphs_probs = p.^sum(all_graphs, 2) .* (1-p).^ (num_edges-sum(all_graphs, 2));

probs_time = cputime - ttt

all_graphs_overlaps = all_graphs * all_graphs'; 
all_graphs_probs = graphs_probs * graphs_probs';

matmult_time = cputime - ttt


x_mean = zeros(1, length(alpha_vec)); x_std = zeros(1, length(alpha_vec));
for i=1:length(alpha_vec)
    now_alpha = alpha_vec(i)
    x_mean(i) = sum(all_graphs_probs(find(all_graphs_overlaps >= alpha_vec(i)*num_edges)));
end

threshes_time = cputime - ttt

% Now loop over all graphs
% % % x_mean = 0; 
% % % for i=1:2^num_edges
% % %     cur_overlap = sum(repmat(all_graphs(i,:),  2^num_edges, 1) .* all_graphs, 2);
% % %     
% % %     x_mean = x_mean + graphs_probs(i) * sum(graphs_probs(find(cur_overlap >= alpha * num_edges)));
% % %     
% % % end

x_mean = x_mean .* factorial(N);



% Now to the challenging part: Compute the variance of X:
gp = graph_perms(N); 

% Prepare a permutation among graphs ! 
big_graph_perms = zeros(2^num_edges, factorial(N)); 

for i=1:factorial(N)
    temp_all_graphs = all_graphs(:,gp(i,:));
    
    big_graph_perms(:,i) = temp_all_graphs * (2.^[0:num_edges-1]') + 1;
     
end


big_perms_time = cputime - ttt

% Measure correlation E X_i X_{\pi}
perms_corrs = zeros(length(alpha_vec)+2, factorial(N));


% Calculate the overlap
nodes_perms = perms([1:N]);
size(nodes_perms)
size(repmat([1:N], factorial(N), 1))
perms_corrs(end-1,:) = sum(nodes_perms == repmat([1:N], factorial(N), 1), 2);

perms_corrs(end,:) = calc_swaps_distance(N); % Calculate distance differently, using swaps


for i=1:length(alpha_vec)

    % Indexes of pairs of graphs with sufficiently good alignment
    overlap_invalid_indexes = find(all_graphs_overlaps < alpha_vec(i)*num_edges);
    overlap_valid_indexes = find(all_graphs_overlaps >= alpha_vec(i)*num_edges);
    overlap_matrix = ones(2^num_edges); overlap_matrix(overlap_invalid_indexes)=0;
    
    
    indexes_i=mod(overlap_valid_indexes-1,2^num_edges)+1; indexes_j=floor((overlap_valid_indexes-1)/(2^num_edges))+0*1;

    % Keep only graphs which are aligned well enough with the first permutation
    all_graphs_probs_valid = all_graphs_probs;
    all_graphs_probs_valid(overlap_invalid_indexes)=0; % Make them zero !!! 
    
    
    x_mean_alter(i) = sum(sum(all_graphs_probs_valid));
    
    all_graphs_at_least_one_good_perm = zeros(2^num_edges);
    
    pp_inv = zeros(1,2^num_edges);
    
        
    
    for(j=1:factorial(N)) % loop over all graph permutations
        
        pp_inv(big_graph_perms(:,j)) = [1:2^num_edges];
        
        cur_perm_good_indexes = ((2^num_edges)*(indexes_j)+(pp_inv(indexes_i)'));
     
   %     cur_perm_good_indexes = find(overlap_matrix(big_graph_perms(:,j),:));
%        cur_perm_good_indexes = find(all_graphs_overlaps(:,big_graph_perms(:,j)) >= alpha_vec(i)*num_edges);
        all_graphs_at_least_one_good_perm(cur_perm_good_indexes)=1;
                
        
        perms_corrs(i,j) = sum(all_graphs_probs_valid(cur_perm_good_indexes));
                
        
        x_std(i) = x_std(i) + perms_corrs(i,j); %%% sum(all_graphs_probs_valid(cur_perm_good_indexes));
        %%%jj=j
    end
            
    p_x_positive_exact(i) = sum(sum(all_graphs_at_least_one_good_perm .* all_graphs_probs));
    
    iii=i
end

x_std = x_std .* factorial(N);


% Now do the 2nd moment method
p_x_positive_2nd_moment = x_mean.^2 ./ x_std;

x_std = sqrt(x_std - x_mean.^2);


x_mean_alter
x_mean ./ factorial(N)










    
    
    
            
    
    

    
    
    



