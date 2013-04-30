% Generate random samples from a BNT (given as a distribution) using mcmc 
% my_mcmc_sample_bnet_counts Generate many random samples from a Bayes net, (which is already given as a distribution)
% of size nsamples, and the number of samples is given by num_iters
%
% Here the sampling is done using monte-carlo markov chain simulations
% sample(i,j) contains the count of the i'th node in the j'th sample vec
% Nodes are sampled in the order given by bnet.order.
% Note : This is a special-purpose sampling function for fast sampling of
% small BINARY discrete BNT !!!! We also assume we have less than 32 nodes,
% so that each sample can fit into one word !!!!!!
% It is not good for the general case.
% Note : Unlike the usual sample_bnet, this function returns an matrix and
% not a cell-array
%
% num_samples is the number of samples to draw from the distribution.
% num_iters is the number of repetitions to do when sampling, so we get
% num_iters different samplings, each of size num_samples.
% We now do it in the new way. Note that now the function returns the
% sufficient statistics, or the sample counts (NOT normalized by the number of samples taken),
% and NOT the entire sample !!!!  This is for now
% good only for small BNETs. Thus we do not need to collect the samples later
function sample_counts = my_mcmc_sample_bnet_counts(init_count, num_mix_steps, bnet_cum_dist, ...
    num_iters, binom_cumsum_vecs, varargin)

% Prepare the ratio's matrix
m = length(bnet_cum_dist); % usually m=2^n

p_ratios_vec = zeros(1,m*m);

for i=1:m
    p_ratios_vec((i-1)*m+1:i*m) = bnet_cum_dist(i) ./ bnet_cum_dist;
end



% Set the initial starting point for the markov chain
sample_count = init_count;

% perform many slack steps to 'mix' the chain
sample_counts = my_mcmc_steps(sample_count,  num_mix_steps, p_ratios_vec);

% Now perform the 'true' steps for sampling
sample_counts = my_mcmc_steps(sample_counts(end,:),  num_iters, p_ratios_vec);



% Here actually perform the steps
function sample_counts = my_mcmc_steps(init_count,  num_steps, p_ratios_vec);



sample_count = init_count;
m=sqrt(length(p_ratios_vec));

% check if this is integer
if(floor(m) < m)
    
    sprintf('Errorrrrrr !!!! NOt integer M !!!!')
    sample_counts = 999999;
    return;
end
    
    
sample_counts = zeros(num_steps, m);


% Change the format for easier randomization
%%%%%%%%%%p_ratios_vec = reshape(p_ratios_mat, 1, m*m);

% For now do it step by step (in the future vectorize)
rand_choose_vec = p_ratios_vec;


rand_vec = rand(1,num_steps);

for iter=1:num_steps

    % Compute whatever we want to compute on this - For now no! We just
    % keep in memory the points !!!


    % Prepare the matrix to draw from
    for i=1:m
        rand_choose_vec((i-1)*m+1:i*m) = p_ratios_vec((i-1)*m+1:i*m) .* ( sample_count(i) ./ (sample_count+1) );
    end

    % Each time compute the cumsum again ???
    rand_choose_cumsum_vec = cumsum(rand_choose_vec);

    % Use binary search routine of Eran Ofek, to find the index of closest value (we need upper ..)
    % avoid 0 index
    ind=max( bin_sear(rand_choose_cumsum_vec,rand_vec(iter)), 1);
    if(rand_choose_cumsum_vec(ind) >= rand_vec(iter))
        C_i = floor((ind-1)/m)+1; C_j = mod(ind-1, m)+1;
    else  % here our rand vec has passed the value so we have to add one
        C_i = floor(ind/m)+1; C_j = mod(ind, m)+1;
    end


    % Do the flip
    sample_count(C_i) = sample_count(C_i)-1;
    sample_count(C_j) = sample_count(C_j)+1;
    
%     sum(sample_count)
%     size(sample_count)
%     size(sample_counts(iter,:))
    
    % Save the count
    sample_counts(iter,:) = sample_count;

end