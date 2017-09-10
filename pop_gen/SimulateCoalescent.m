% Simple simulation of coalescent process with no recombination
% Input:
% n - number of individuals
% theta - scaled mutation rate
% M - number of simulations to perform
%
% Output:
% Trees - structure with information on simulated trees
% Mut - structure with infromation on simulated mutations
%
function [Mut, Trees] = SimulateCoalescent(n, theta, M)

lambda_vec = 0.5 .* (2:n) .* (1:n-1); % set rate at each stage

T_mat = zeros(M, n-1);   % set coalescent times
for i=2:n
    T_mat(:,i-1) = exprnd(1./lambda_vec(i-1), M, 1); % Matlab exprnd requires 1/lambda
end

% Set mutations
T_MRCA = sum(T_mat, 2);
T_Total = sum(repmat(2:n, M, 1) .* T_mat, 2);
T_cum_mat = [zeros(M, 1) cumsum(T_mat(:,end:-1:1),2)]; T_cum_mat = T_cum_mat(:,end:-1:1);
Trees = [];
Trees.T_MRCA_mean = mean(T_MRCA);
Trees.T_Total_mean = mean(T_Total);

Trees.Topology = zeros(M, 2*(n-1)); % New: keep also tree identity
for i=1:(n-1) % loop on coalescent events
    Trees.Topology(:, 2*i-1) = ceil(rand(M,1) .* (n+1-i)); % choose first branch to coaelesce
    Trees.Topology(:, 2*i) = ceil(rand(M,1) .* (n-i));  % choose second branch to coaelesce
    tmp_inds = Trees.Topology(:, 2*i) >= Trees.Topology(:, 2*i-1);
    Trees.Topology(tmp_inds, 2*i) = Trees.Topology(tmp_inds, 2*i) + 1;
end


Mut = [];
Mut.num = poissrnd(T_Total .* theta ./ 2); % set # of mutations (sites) in each simulation
Mut.levels = cell(M,1); Mut.ages = cell(M,1); Mut.branches = cell(M,1); Mut.counts = cell(M,1);
Mut.n_alleles = zeros(M,1);  Mut.n_alleles2 = zeros(M,1);
for i=1:M % do non-vector (might improve in future
    if(mod(i, 1000) == 0)
        run_i = i % plot to see advances in running 
    end
    Mut.levels{i} = weighted_rand(T_mat(i,:) .* (2:n), Mut.num(i))+1;
    Mut.ages{i} = rand(1, Mut.num(i)) .* T_mat(i, Mut.levels{i}-1) + T_cum_mat(i, Mut.levels{i}); % set time within the level plus time from previous levels
    Mut.branches{i} = ceil(rand(1, Mut.num(i)) .* Mut.levels{i}); % Set mutation branches (index)
    
    G = genotype_mat_from_indices(Trees.Topology(i,:), Mut.levels{i}, Mut.branches{i});
    Mut.n_alleles2(i) = length(unique(Mut.levels{i}))+1; % count mutations on same branch only once !!
    Mut.n_alleles(i) = size(unique(G, 'rows'), 1); % Count # alleles
    
    TmpTree = tree_from_pairwise_indices(Trees.Topology(i,:)); % get tree
    % get allele freq.
    Mut.counts{i} = zeros(size(Mut.levels{i}));
    for j=1:length(Mut.levels{i})
        Mut.counts{i}(j) = length(TmpTree{n+1-Mut.levels{i}(j)}{Mut.branches{i}(j)});
    end
end

% Compute summary statistics
Mut.age_mean = mean(cell2vec(Mut.ages));
Mut.age_var = var(cell2vec(Mut.ages));
Mut.levels_hist = hist(cell2vec(Mut.levels), 2:n);
Mut.levels_freq = Mut.levels_hist ./ sum(Mut.levels_hist);

