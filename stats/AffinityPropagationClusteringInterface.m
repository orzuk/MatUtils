% Cluster vectors using Frey's affinity-propagation algorithm. 
% Before clustering, we need first to compute pair-wise similarity values
% This function merely serves as an interface: 
% First - compute the knn for each point using either ANN or brute-force
% Second - call the apcluster function to do the clustering work for you
% The size of problem which currently can be handled is on the order of
% thousands to tens of thousands 
%
% Input:
% V - the set of vectors that we want to cluster
% counts - the number of times each vector appears (this enables getting many vectors together)
% metric - which metric to use
% knn - number of nearest neighbors to keep
% alpha - the 'clustering parameter'. Higher alpha means more clusters. Zero means min similarity and 0.5 means median similarity
% outfile - where to save the output (put empty if you don't want to save output) 
%
% Output: 
% num_clusters - number of obtained clusters
% cluster_sizes - size of each cluster
% cluster_idx - the indices of each point, representing to which cluster does it belong
% netsim - similarity scores achieved by the clustering 
% cluster_centers - means of clusters (rather than examplars)
%
function [num_clusters cluster_sizes cluster_idx netsim cluster_centers] = ...
    AffinityPropagationClusteringInterface(V, counts, metric, knn, alpha, outfile, varargin)

if(iscell(V)) % for cell-arrays we need the correct orientation 
    if(size(V,1) == 1)
        V = V';
    end
end

n = size(V, 2); % dimension of vectors
m = size(V, 1); % number of vectors we want to cluster

if( (~exist('knn', 'var')) && isempty(knn) )
    knn = m;
end
knn = min(m, knn); % maximum number of neighbors

ttt_knn = cputime;
% First calculate the knn by brute-force. this should be later replaced by a better algorithm
[KNN_sim KNN_inds] = KNNBruteForce(V, [], metric, knn); 
if(exist('outfile', 'var') && (~isempty(outfile)))
    my_mkdir(dir_from_file_name(outfile)); % make sure directory exist
    save(outfile, 'V', 'KNN_sim', 'KNN_inds'); % save knn results and also the vectors 
end
ttt_knn = cputime - ttt_knn % display just the knn computation time
 
% KNN_dists = sparse(KNN_dists); 
% [I J S] = find(KNN_dists); % generate the distance matrix S
% [I_inds J_inds S_inds] = find(KNN_inds); 

ttt_apcluster = cputime;
p = zeros(m,1) + my_quantile(KNN_sim(:), alpha); % preferences value - determines roughly the number of clusters

if(~isempty(counts)) % adjust for multiplicity
    %        KNN_sim = KNN_sim .* repmat(counts, 1, m) .* repmat(counts', m, 1); % sim(i,j) * c(i) * c(j)
    p = p + (counts-1) .* KNN_sim(:,1);  % p(i) + (c(i)-1)*sim(i,i)
    KNN_sim = KNN_sim .* repmat(counts, 1, knn); % sim(i,j) * c(i)
    %       for j=1:knn
    %           KNN_sim(:,j) = KNN_sim(:,j) .* counts(KNN_inds(:,j)); % multiplying also by c(j) - NOT NEEDED !!! (only one examplar)
    %       end
 end
S = [ repmat(1:m, 1, knn)'  reshape(KNN_inds, 1, m*knn)'  reshape(KNN_sim, 1, m*knn)' ]; % I inds, J inds, distances

% % M=N*N-N; s=zeros(M,3); % Make ALL N^2-N similarities
% % j=1;
% % for i=1:N
% %     for k=[1:i-1,i+1:N]
% %         s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);
% %         j=j+1;
% %     end;
% % end;
% % p=median(s(:,3)); % Set preference to median similarity

tic; [cluster_idx,netsim,dpsim,expref] = apcluster(S,p); toc % % apclustersparse (doesn't work well)  apclustersparsemex(S,p); Now call the apcluster function
cluster_idx = max(cluster_idx, 1); % handles a bizarre case in which we get 0 index 
num_clusters = length(unique(cluster_idx));
cluster_sizes = hist(double(cluster_idx), unique(double(cluster_idx)));
tic; cluster_centers = get_cluster_centers(V', cluster_idx); toc


ttt_apcluster = cputime - ttt_apcluster


cluster_again_flag = 0; % flag saying if to cluster again ... 
if(cluster_again_flag == 1) % cluster another level ...
    u = unique(cluster_idx);
    
    V = V(u, :);
    [KNN_sim KNN_inds] = KNNBruteForce(V, [], metric, knn); 
    p = zeros(length(u),1) + my_quantile(KNN_sim(:), alpha); % preferences value - determines roughly the number of clusters
    counts = cluster_sizes';
    p = p + (counts-1) .* KNN_sim(:,1);  % p(i) + (c(i)-1)*sim(i,i)
    KNN_sim = KNN_sim .* repmat(counts, 1, knn); % sim(i,j) * c(i)

    S = [ repmat(1:length(u), 1, knn)'  reshape(KNN_inds, 1, length(u)*knn)'  reshape(KNN_sim, 1, length(u)*knn)' ]; % I inds, J inds, distances
    tic; [cluster_idx2,netsim,dpsim,expref] = apclustersparsemex(S,p); toc % Now call the apcluster function
    num_clusters2 = length(unique(cluster_idx2));
    cluster_sizes2 = hist(double(cluster_idx2), unique(double(cluster_idx2)));
    true_cluster_sizes2 = zeros(1,length(cluster_sizes2));
    u2 = unique(cluster_idx2); 
    for j=1:length(u2)
        u2_inv_perm(u2(j)) = j;
    end
    for j=1:length(u)
        true_cluster_sizes2(u2_inv_perm(cluster_idx2(j))) = true_cluster_sizes2(u2_inv_perm(cluster_idx2(j))) + counts(j);
    end
    [mm mmm] = max(true_cluster_sizes2)
    cluster_centers2 = get_cluster_centers(V', cluster_idx2);
end

if(exist('outfile', 'var') && (~isempty(outfile)))
    save(outfile, 'cluster_idx', 'netsim', 'dpsim', 'expref', 'num_clusters', 'cluster_sizes', '-append'); % save results (keep the knn stuff)
end
