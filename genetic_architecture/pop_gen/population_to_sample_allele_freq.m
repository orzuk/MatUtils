% Take a random sample from a population and compute allele frequencies in sample
% Input:
% f_vec - frequency of each allele
% N - population size - total number of CHROMOSOMES !!!!
% n_vec - number of individual CHROMOSOMES!!! in sample for each allele (usually all the same)
% return_counts - flag saying to return counts (default) or frequencies 
%
% Output:
% k_vec - number of allele carriers for each allele (drawn at random)
%
function k_vec = population_to_sample_allele_freq(f_vec, N, n_vec, return_counts)

if(~exist('return_counts', 'var') || isempty(return_counts))
    return_counts = 1;
end
f_vec = vec2column(f_vec);
num_alleles = length(f_vec);
poiss_flag=1; % approximation 
if(isscalar(n_vec))
    if(n_vec < 0.1*N)  % If n << N we use Poisson approximation 
        poiss_flag=1;
    end
    n_vec = repmat(n_vec, num_alleles, 1);
end

if((max(f_vec) <= 1.001) && ~isempty(find(f_vec(f_vec>0) < 1, 1))) % transfer alleles frequencies to counts
    f_vec = round(f_vec .* N);
end

if(poiss_flag) % might save time 
%    poiss_inds = find(max(n_vec, round(f_vec)) < 0.1*N); hyg_inds = setdiff(1:num_alleles, poiss_inds);
%    k_vec = zeros(num_alleles, 1); 
%    k_vec(poiss_inds) = poissrnd( vec2column(round(f_vec(poiss_inds))) .* n_vec(poiss_inds) ./ N); % much faster: use poisson approximation
%    k_vec(hyg_inds) = hygernd(N, vec2column(round(f_vec(hyg_inds))), n_vec(hyg_inds)); % Sample without replacement. Can sample together alleles with same frequencies? probably NOT!
%    save_f_vec = f_vec; 
    big_inds = find(f_vec > N/2); small_inds = setdiff(1:num_alleles, big_inds); k_vec = zeros(num_alleles, 1); f_vec(big_inds) = N-f_vec(big_inds); 
    mu =  f_vec .* n_vec ./ N; 
    norm_inds = find(mu > 50); % find(f_vec > 1000); %  & ((N-f_vec)>1000)); 
    hyg_inds = setdiff(1:num_alleles, norm_inds); 
    k_vec(hyg_inds) = hygernd(N, vec2column(round(f_vec(hyg_inds))), n_vec(hyg_inds));
    k_vec(norm_inds) = round( normrnd( f_vec(norm_inds) .* n_vec(norm_inds) ./ N, ...
        sqrt(n_vec(norm_inds) .* (N-n_vec(norm_inds)) .* f_vec(norm_inds) .* (N-f_vec(norm_inds)) ./ (N*N*(N-1)) )) ); % Use Gaussian approximation 
    k_vec(norm_inds) = min(max(0, k_vec(norm_inds)), n_vec(norm_inds)); 
    k_vec(big_inds) = n_vec(big_inds) - k_vec(big_inds); 
%    f_vec = save_f_vec;
    
%     big_inds = intersect(big_inds, hyg_inds); small_inds = intersect(small_inds, hyg_inds); 
%     k_vec(big_inds) = n_vec(big_inds) - hygernd(N, N-vec2column(round(f_vec(big_inds))), n_vec(big_inds));
%     k_vec(small_inds) = hygernd(N, vec2column(round(f_vec(small_inds))), n_vec(small_inds));
%     mu =  f_vec(norm_inds) .* n_vec(norm_inds) ./ N; 
%     sigma = sqrt(n_vec(norm_inds) .* (N-n_vec(norm_inds)) .* f_vec(norm_inds) .* (N-f_vec(norm_inds)) ./ (N*N*(N-1)));
%     k_vec(norm_inds) = round( normrnd( mu, sigma )); % Use Gaussian approximation 
%     poiss_inds = norm_inds(mu < 50); norm_inds = setdiff(norm_inds, poiss_inds);  
%     k_vec(poiss_inds) = poissrnd( vec2column(round(f_vec(poiss_inds))) .* n_vec(poiss_inds) ./ N); % much faster: use poisson approximation
else  % exact, can be slow
    k_vec = hygernd(N, vec2column(round(f_vec)), n_vec);
end
if(~return_counts) % return frequencies rather than counts (counts is default)
    k_vec = k_vec ./ n_vec;
end


