% Generate a graph according to some model.
% Then sample from it and see how many things change
% due to the sampling process

function [f_mat corr_mat] = check_ranking_stability(N, alpha, iters, p_vec1, p_vec2, samp_type, g_t, param)
SAMP_EDGES = 0; SAMP_NODES = 1;  % sampling flags

TOL = 0.000000001;

% First randomize a graph
G=generate_random_graph( g_t, N, param); deg_vec = full(sum(G));

[deg_val deg_ind] = sort(-deg_vec); % - indicates we want the high degrees

%rp = randperm(N); deg_val = deg_val(rp); deg_ind = deg_ind(rp);


deg_inv_ind = inv_perm(deg_ind, N); % Inverse the permutation

deg_val = deg_val + TOL*rand(1,length(deg_val));

top_vars = round(alpha * N);

% Permute the ties to be fair
 deg_x_alpha = -deg_val(top_vars);
% deg_x_alpha_inds = find(deg_val == deg_x_alpha);
% rp = randperm(length(deg_x_alpha_inds));
% deg_ind(deg_x_alpha_inds) = deg_ind(deg_x_alpha_inds(rp)); 


f_mat = zeros(length(p_vec1), iters);

% New! Calculate also the correlation between the permutations - a
% different measure which is independent of alpha
corr_mat = zeros(length(p_vec1), iters);

% p_ind=1;

edge_dense = 0;

for p_ind=1:length(p_vec1)
    p1 = p_vec1(p_ind); p2 = p_vec2(p_ind);
    for i=1:iters
        Gs = sample_graph(G, samp_type, p1, p2); deg_vec_s = full(sum(Gs)) + TOL*rand(1,length(deg_val));

        
        edge_dense = edge_dense + full(sum(sum(Gs)));
        
        % Now compare the overlap to the original
        [deg_val_s deg_ind_s] = sort(-deg_vec_s); % - indicates we want the high degrees
        
        deg_inv_ind_s = inv_perm(deg_ind_s, N); % Inverse the permutation
        
        % Note: Here we need to cleverly deal with ties!
        f_mat(p_ind, i) = length(intersect(deg_ind(1:top_vars),deg_ind_s(1:top_vars)));

        
%        aaa = size(corr(deg_ind, deg_ind_s))
        corr_mat(p_ind,i) = corr(deg_inv_ind', deg_inv_ind_s');
        
        
    end
%     p_ind=p_ind+1;
end


edge_dense = edge_dense / (iters*length(p_vec1));


edge_bias = edge_dense - full(sum(sum(G)))

% Normalize to get overlap between zero and one
f_mat = f_mat ./ top_vars; 


% % Do some figures
% figure; subplot(2,2,1); hold on; errorbar(p_vec, mean(f_mat,2), std(f_mat,[],2));
% title(['overlap in top degrees \alpha=' num2str(alpha) ' N = ' num2str(N)]); xlabel('p'); ylabel('f');
% subplot(2,2,2); plot(p_vec, mean(f_mat,2));
% 
% % Plot the loglog degreee distribution:

% % subplot(2,2,3); hold on;  
% r=full(sum(G));  h = hist(r, max(r)-min(r)+1);
% figure;  xlabel('k'); ylabel('P(k)'); loglog(min(r):max(r), h./N, '*');   
% %loglog( deg_x_alpha  , 0, '*r'); % This is the cutoff degree
% 
% % Plot histogram of the degree distribution
% deg_un = size(unique(r),2)
% figure; hold on;
% hist(r, deg_un); title('G Degree dist. ');


function ip = inv_perm(p, N)
for i=1:N
    ip(p(i))=i;
end


