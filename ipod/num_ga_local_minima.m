% Compute number of local minima of alignments between two graphs.
% The function gets as input two graphs, G1 and G2,
% and computes the number of local minima, i.e. locally optimal
% alignments between the two graphs.
% 
function [global_minima local_minima ground_state_energy] = num_ga_local_minima(G1,G2)


N = size(G1,1);

if(N > 10) % Too big N
    sprintf('Error! Cannot exhaustively search aligments for N>10\n')
    global_minima=[];local_minima = [];
end

all_ga_scores = find_all_graph_alignments(G1,G2); ground_state_energy = max(all_ga_scores);
global_minima = sum(all_ga_scores == ground_state_energy);
all_perms = perms(1:N);
all_perms = sortrows(all_perms-1);
local_minima_vec = ones(1,factorial(N));

% Generate all transpositions
all_trans = repmat(1:N,N*(N-1)/2,1);
i=1;
for j=1:N
    for k=j+1:N
        all_trans(i,j)=k; all_trans(i,k)=j; i=i+1;
    end
end

% all_trans=all_trans-1
% size(all_trans)


% Go over all permutations
for i=1:factorial(N)
    if(local_minima_vec(i))
        % Take the current permutation and generate all of its swaps
        swap_perm = all_perms(i,:);
        swap_perm = swap_perm(all_trans);
        swap_index = perm_to_index(swap_perm);
        if(max(all_ga_scores(swap_index)) > all_ga_scores(i))
            local_minima_vec(i) = 0;
        end
        local_minima_vec( swap_index( find( all_ga_scores(swap_index)  < all_ga_scores(i) ) ) ) = 0;
        if(mod(i,1000)==0)
            doing_i = i
        end
    end
end

% local_minimas_are = find(local_minima_vec)
% minima_scores = all_ga_scores(local_minimas_are)
local_minima = sum(local_minima_vec);


