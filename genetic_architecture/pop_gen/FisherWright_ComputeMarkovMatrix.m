% Compute the transition matrix for one step of FisherWright model
% Input:
% N - population size
% s - selection coefficient
% prob_compute_method - how to compute transition probabilities. 'exact' - binomial probabilities, 'approximate' - use Poisson or Normal approximation)
% compute_matrix - flag saying if to compute the entire matrix, or just advance to compute next p
% p_cur - (optional) current allele frequency probabilities
%
% Output:
% M - transition matrix
% p_next - probability distribution at the next generation
%
function [M, p_next] = FisherWright_ComputeMarkovMatrix(N, s, prob_compute_method, compute_matrix, compute_prob_vector,  p_cur)

if(~exist('s', 'var') || isempty(s)) % default is neutral
    s=0;
end
if(~exist('prob_compute_method', 'var') || isempty(prob_compute_method))
    prob_compute_method = 'approximate'; % 'exact'; % 'approximate'; %'exact'; % s'exact'; % ''; % 'exact';
end
if(isscalar(N)) % set current and next generation population size
    N = [N N];
end
if(~exist('compute_matrix', 'var') || isempty(compute_matrix))
    compute_matrix = 1; % compute transition matrix
end
if(~exist('compute_prob_vector', 'var') || isempty(compute_prob_vector))
    compute_prob_vector = 0; % compute next generation probability vector for frequencies 
end
x_vec = (0:2*N(1)) ./ (2*N(1)); % set current allele frequency vector

new_x_mean_vec = x_vec .* (1+s) ./ (1 + x_vec.*s); % first get the mean vector of new #offspring
num_sigmas = 5; % keep only these many st.d. (to speed-up computation)
std_vec = sqrt(2*N(1) .* new_x_mean_vec .* (1-new_x_mean_vec)); % binomial standard deviation

if(compute_prob_vector)
    p_next = zeros(2*N(2)+1, 1); % vector of probabilities in next generation
end
if(compute_matrix)
    M = zeros(2*N(1)+1, 2*N(2)+1); % allow change in population size
end

left_range_vec = max(2, round(new_x_mean_vec.*2*N(2) - num_sigmas .* std_vec)); % NEW! start from 2, not 1 !!!! % calculate range to compute binomial distribution
right_range_vec = 1+min(2*N(2), round(new_x_mean_vec.*2*N(2) + num_sigmas .* std_vec));
if(compute_prob_vector)
    non_negligile_x_inds = find(p_cur > 10^(-9)); % take only states with non-negligible probability
end
if(compute_matrix)
    run_x_inds = 1:(2*N(1)+1); 
else
    run_x_inds = non_negligile_x_inds;
end
for k=vec2row(run_x_inds) % 1:2*N(1)+1 % -1 % loop on current allele frequency. This is the heaviest part
    cur_range_vec = left_range_vec(k):right_range_vec(k);
    
    % Exact binomial computation
    switch prob_compute_method
        case 'exact'
            if(compute_prob_vector)
                p_next(cur_range_vec) = ...
                    p_next(cur_range_vec) + p_cur(k)  .* ...
                    binopdf(cur_range_vec-1, 2*N(2), new_x_mean_vec(k)); % why take cur_range_vec - 1 ???
            end
            if(compute_matrix)
                M(k,:) = binopdf(0:2*N(2), 2*N(2), new_x_mean_vec(k)); % Why column here? % -1
            end
        case {'approx', 'approximate'}
            if(new_x_mean_vec(k) * 2*N(2) < 50)  % use poisson approximation  to binomial on left tail
                if(compute_prob_vector)
                    p_next(cur_range_vec) = ...
                        p_next(cur_range_vec) + p_cur(k)  .* ...
                        poisspdf(cur_range_vec-1, 2*N(2) * new_x_mean_vec(k));  % Replace with poisson approximation for larger values!!!
                end
                if(compute_matrix)
                    M(k,:) = poisspdf(0:2*N(2), 2*N(2) * new_x_mean_vec(k)); % Poisson approximation
                end
            else % for larger values use poisson or Gaussian approximations
                if((1-new_x_mean_vec(k)) .* 2*N(2) < 50) %  use poisson approximation to binomial on right tail,
                    if(compute_prob_vector)
                        p_next(cur_range_vec) = ...
                            p_next(cur_range_vec) + p_cur(k)  .* ...
                            poisspdf(2*N(2) - (cur_range_vec-1), 2*N(2) * (1-new_x_mean_vec(k)));  % This part is slowest. We replaced with poisson !!!
                    end
                    if(compute_matrix)
                        M(k,:) = poisspdf(2*N(2):-1:0, 2*N(2) * (1-new_x_mean_vec(k))); % Poisson approximation
                    end
                else % use Gaussian approximation to binomial in the middle of distribution
                    if(compute_prob_vector)
                        p_next(cur_range_vec) = ...
                            p_next(cur_range_vec) + p_cur(k)  .* ...
                            normpdf_vec( round(( 6+(cur_range_vec-1 - 2*N(2) * new_x_mean_vec(k)) ./ std_vec(k) ) .* 10^5) ) ./ std_vec(k);
                    end
                    %                 normpdf(cur_range_vec-1, 2*N(2) * new_x_mean_vec(k), ...
                    %                 sqrt(2*N(2) * new_x_mean_vec(k) * (1-new_x_mean_vec(k))) );  % This part is slowest. We replace with poisson !!!
                    if(compute_matrix)
                        M(k,:) = normpdf(0:2*N(2), 2*N(2)*new_x_mean_vec(k), ...
                            sqrt(2*N(2)*new_x_mean_vec(k)* (1-new_x_mean_vec(k))) ); % Why column here? % -1
                    end
                end
            end % if: what probability propagating approximation to use
    end % switch approximate
    
end % loop on k


