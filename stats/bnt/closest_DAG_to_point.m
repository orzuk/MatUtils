% Find the DAG (model) which is closest in relative entropy to
% a given distribution (point). We allow the maximum in-degree
% to be max_k. We currently check only DAGs with the  'natural' orientation.
% Checking all orientations is exponential.
function [G KL] = closest_DAG_to_point(P, max_k)

KL_MI = 0;
n = log2(size(P,2)); % Get the dimenstion
G = zeros(n);
max_k = min(max_k, n-1); % We cannot have more than n-1 edges

% For the first k-1, we simply generate the click
for i=1:max_k+1
    G(1:i-1,i) = 1;
end

% Loop over variables and find the set of parents getting us
% closest to P for each variable
for i=max_k+2:n
    I_is = i;
    candidate_parents_vec = nchoosek(1:i-1,max_k);    % generated all k-tuples
    % Loop over all k-tuples
    best_MI = 99999999999999; best_parents = [];
    for j=1:nchoosek(i-1,max_k)
        cur_MI = conditional_mutual_information(P, i, setdiff(1:i-1, candidate_parents_vec(j,:)), candidate_parents_vec(j,:));
%         if(cur_MI < best_MI)
%             best_MI = cur_MI;
%             best_parents = candidate_parents_vec(j,:);
%         end
        best_MI = min(cur_MI, best_MI);

    end

%     G(best_parents,i) = 1; % update the graph G
    KL_MI = KL_MI + best_MI;
end

% Get the relative entropy distance
% KL_MI_is = KL_MI
KL = KL_MI;
% KL2 = relative_entropy_to_bnet(P', G);
% diff_is = KL2 - KL;
% ratio_is = KL2 / KL;
