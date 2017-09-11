% Compute statistics for several architectures
% For now work only for binary architectures
%
function [v_pairwise v_pairwise_explained lods_ratio_pairwise mu_pairwise] = ...
    compute_architecture_statistics_pairwise(architecture_str, ...
    V, f_vec, params_struct, p_x_vec, p_x_times_z_vec, x_ind_mat, iters, ...
    v_marginal, lods_ratio_marginal, ...
    compute_method_flag, compute_pairwise_flag) % compute several moments and other stuff for architecture

AssignGeneralConstants;
AssignStatsConstants;
N = length(f_vec);
% v_pairwise = zeros(N); lods_ratio_pairwise = zeros(N,N,4);  

if( ~exist('x_vec', 'var') || isempty(x_vec)) % initilize x_vec if needed
    switch compute_method_flag
        case 'sampling' % sample a few x's
            %    x_vec = double(rand(iters,N) < 0.5);
            x_vec = double(rand(iters,N) < repmat(f_vec, iters, 1)); % note: used here the same iters for two things !!!
            %        iters = min(100, iters); % lower sampling from multiple x's !!!
        case 'enumerate' % enumerate all 2^N possible vectors
            x_vec =  dec2bin(0:2^N-1) - '0'; % take all 2^N possibilities
            
        case 'analytic' % analyic computation. No need to enumerate
    end
end


arch_type = architecture_to_type(architecture_str);
switch compute_method_flag % how to compute
    case 'enumerate'
%         [z z_clean] = ...
%             genetic_architecture(x_vec, architecture_str, params_struct, iters); % generate outputs
        %        V = sum(p_x_vec .* z_clean.^2) / iters - sum(p_x_vec .* z_clean)^2 / iters^2 + params_struct.z_std^2; % variance of trait
        iters = 1; v_pairwise = zeros(N);
        
        if(compute_pairwise_flag)
            mu_pairwise = zeros(N,N,4); 
            ttt_pairwise = cputime;
            p_ij = zeros(2);
            for i=1:N
                for j=i+1:N
                    ind_vec = zeros(N,1); ind_vec(i) = 1; ind_vec(j) = 1;
                    v_pairwise(i,j) = ...
                        conditional_var(x_vec, p_x_vec, p_x_times_z_vec, ind_vec, iters); % V(z|x_i,x_j)
                    for x_i=0:1
                        %            x_i_ind_vec = find(x_vec(:,i) == x_i);
                        for x_j = 0:1
                            %                x_j_ind_vec = find(x_vec(:,j) == x_j);
                            %                cur_ind_vec = intersect(x_ind_vec{i,x_i+1}, x_ind_vec{j,x_j+1});
                            p_ij(x_i+1,x_j+1) = ...
                                sum( p_x_times_z_vec(x_ind_mat{i,j,x_i+1,x_j+1}) ) ...
                                / sum(p_x_vec(x_ind_mat{i,j,x_i+1,x_j+1}));  % Pr(z=1 | x_i,x_j)
                        end
                    end
                    mu_pairwise(i,j,:) = mat_into_vec(p_ij ./ p_ij(1,1)); % normalize by the baseline risk of (0,0)
                end
            end
%             for i=1:4% Now correct for the individual log-ratios
%                 lods_ratio_pairwise(:,:,i) = lods_ratio_pairwise(:,:,i) ./ ...
%                     (repmat(lods_ratio_marginal, 1, N) .* repmat(lods_ratio_marginal, 1, N)');
%             end
%             lods_ratio_pairwise = lods_ratio_pairwise(:,:,4) ./ lods_ratio_pairwise(:,:,1);
            ttt_pairwise = cputime - ttt_pairwise;
        end
        
    case 'sampling'
    case 'analytic' % works for certain architectures        
        iters = length(params_struct.z_std); 
        if(compute_pairwise_flag) % compute pairwise
            mu_pairwise = zeros(N,N,iters,4); % Pr(z=1|x_i,x_j)
            v_pairwise = zeros(N,N,iters,4); % Var(z=1|x_i,x_j)
            v_environment_pairwise = zeros(N,N,iters,4);
            for i=1:N % loop on individual variables and use recursion to compute their
                for j=i+1:N
                    f_vec_pairwise = f_vec;
                    for x_i = 0:1
                        for x_j = 0:1
                            four_ind = x_i*2 + x_j + 1;
                            f_vec_pairwise(i) = x_i; f_vec_pairwise(j) = x_j;
                            [mu_pairwise(i,j,:,four_ind) v_pairwise(i,j,:,four_ind) ...
                                v_environment_pairwise(i,j,:,four_ind)] = ...
                                compute_architecture_statistics(architecture_str, ...
                                f_vec_pairwise, params_struct, p_x_vec, p_x_times_z_vec, ...
                                iters, compute_method_flag); % compute several moments and other stuff for architecture
                        end
                    end
                end
            end
            f_mat = reshape(repmat(f_vec' * f_vec, iters, 1), N, N, iters); % matrix of pairwise probabilities
            f_repmat = reshape(repmat(f_vec, N, iters), N, N, iters);
            f_repmat_transpose = reshape(repmat(f_vec', iters, N)', N, N, iters);
            v_pairwise = v_pairwise(:,:,:,1) .* (1+f_mat - f_repmat - f_repmat_transpose) + ...
                v_pairwise(:,:,:,2) .* (f_repmat_transpose - f_mat) + ...
                v_pairwise(:,:,:,3) .* (f_repmat - f_mat) + ...
                v_pairwise(:,:,:,4) .* f_mat;            
        end
end

v_pairwise_explained = zeros(N,N,iters);
for i=1:iters % here loop .. (not vectorized)
    v_pairwise(:,:,i) = v_pairwise(:,:,i) + v_pairwise(:,:,i)'; % make symmetric. This is Var(z | x_i, x_j)
    v_pairwise_explained(:,:,i) = repmat(v_marginal(:,i), 1, N) + ...
        repmat(v_marginal(:,i), 1, N)' - ...
        V(i) - v_pairwise(:,:,i); % This is Var(z | x_i) + Var(z | x_j) - Var(z) - Var(z|x_i,x_j)
    v_pairwise_explained(:,:,i) = v_pairwise_explained(:,:,i) - ...
        diag(diag(v_pairwise_explained(:,:,i))); % no self interactions !
end

for i=1:N % make symmetric
    for j=i+1:N
        mu_pairwise(j,i,:,:) = mu_pairwise(i,j,:,:);
    end
end
%     mu_pairwise(:,:,i) = mu_pairwise(:,:,i) + mu_pairwise(:,:,i)';
% end
lods_ratio_pairwise = mu_pairwise(:,:,:,4) ./ mu_pairwise(:,:,:,1); %P(z=1|x_i=1,x_j=1) / P(z=1|x_i=0,x_j=0)
for i=1:iters % here loop .. (not vectorized)
    lods_ratio_pairwise(:,:,i) = lods_ratio_pairwise(:,:,i) ./ ...
        (repmat(lods_ratio_marginal(:,i), 1, N) .* ...
        repmat(lods_ratio_marginal(:,i), 1, N)'); % Correct for individual log-ratios
%     lods_ratio_pairwise(:,:,i) = lods_ratio_pairwise(:,:,i) + ...
%         lods_ratio_pairwise(:,:,i)'; % make symmetric.
end

for i=1:N % put marginals on the diagonal 
    lods_ratio_pairwise(i,i,:) = lods_ratio_marginal(i,:); 
end
