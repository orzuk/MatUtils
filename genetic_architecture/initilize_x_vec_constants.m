% Initilize constants such as x vector and its probabilities
%
% Input:
% N - number of loci
% dispose_prob_frac - ignore a certain part of the x's probability space
% f - marginal frequencies
% compute_method_flag - how to compute (sampling, analytic etc.)
% iters - number of different genotypes to sample
%
% Output:
% x_vec - vector of genotypes (binary)
% p_x_vec - vector of probabilities Pr(x = i)
% x_ind_vec - vector of indices (???)
% x_ind_mat - matrix of indices 
%
function [x_vec p_x_vec x_ind_vec x_ind_mat] = ...
    initilize_x_vec_constants(N, dispose_prob_frac, f, ...
    compute_method_flag, iters)

if(~exist('compute_method_flag', 'var') || isempty(compute_method_flag))
    compute_method_flag = 'enumerate';
end
switch compute_method_flag
    case 'enumerate'
        x_vec =  dec2bin(0:2^N-1) - '0'; % take all 2^N possibilities
    case 'sampling' % take a random set of x_vec with probabilities based on f
        x_vec = rand(iters, N) < repmat(f, iters, 1);
end
if(dispose_prob_frac > 0) % throw away many of the low probability x's
    [p_sorted sort_perm] = sort(p_x_vec);
    p_sorted = cumsum(p_sorted);
    last_ind = find(p_sorted > dispose_prob_frac, 1);
    x_vec = x_vec(sort_perm(last_ind:end),:);
    %    p_x_vec = p_x_vec(sort_perm(last_ind:end));
end
p_x_vec = x_to_prob_vec(x_vec, f); % Compute probabilities for all x possibilities

if(nargout > 2) % compute only as needed % M = size(x_vec,1);
    x_ind_vec = cell(N,2); % indices of all variables
    for i=1:N
        for j=1:2
            x_ind_vec{i,j} = find(x_vec(:,i) == j-1);
        end
    end
    x_ind_mat = cell(N,N,2,2);
    for i=1:N
        for j=1:N
            for x_i=0:1
                for x_j=0:1
                    x_ind_mat{i,j,x_i+1,x_j+1} = ...
                        intersect(x_ind_vec{i,x_i+1}, x_ind_vec{j,x_j+1});
                end
            end
        end
    end
end


