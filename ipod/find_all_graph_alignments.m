% Compute all N! possible alignments between two graphs. 
% The function gets as input two graphs, G1 and G2,
% And computes for all the N! permutations, the score
% of aligning these two graphs according to these permutations.
% The function uses exhaustive search therefore it is good only for very
% small graphs!!!!
% 
% Input: 
% G1 - first graph
% G2 - second graph
% 
% Output: 
% all_ga_scores - scores of all possible alignments
% 
function all_ga_scores = find_all_graph_alignments(G1,G2)

N = size(G1,1);

if(N > 10) % Too big N
    sprintf('Error! Cannot exhaustively search aligments for N>10\n')
    all_ga_scores=[]; 
end

all_perms = perms(1:N); 

all_ga_scores = zeros(1,factorial(N));

for i=1:factorial(N)
    all_ga_scores(i) = sum(sum(G1 .* G2(all_perms(i,:),all_perms(i,:))));
end


