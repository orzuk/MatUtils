% Prepare BNT to be sampled. Run only once as pre-processing. 
% my_sample_bnet_prepare Prepares the bnet to be sampled. This is done once, 
% since preperation is quite heavy. Afterwards, we can call my_sample_bnet_generate_sample many times, 
% and get our samples quickly. 
% SAMPLE = SAMPLE_BNET(BNET, ...)
%
% sample(i,j) contains the value of the i'th node in the j'th sample
% Nodes are sampled in the order given by bnet.order.
% Note : This is a special-purpose sampling function for fast sampling of
% small BINARY discrete BNT !!!! We also assume we have less than 32 nodes,
% so that each sample can fit into one word !!!!!!
% It is not good for the general case.
% Note : Unlike the usual sample_bnet, this function returns an matrix and
% not a cell-array
%
function [bnet_cum_dist, binom_probs,binom_cumsum_vecs] = ...
    my_sample_bnet_prepare(bnet, num_samples, chunk_size, varargin)

% Now start sampling. For now we do it very simply : 
% Transfer the BNT to probability distribution, and then sample from this
% distribution
bnet_dist = bnet;  % We now get distributions as inputs, not just bnets ... 
bnet_cum_dist = cumsum(bnet_dist); % We return the cumulative distribution


% We now do it differently (for small number of nodes). Generate a binomial
% vector for each probability
bnet_reverse_cum_dist = cumsum(bnet_dist(end:-1:1));
bnet_reverse_cum_dist = bnet_reverse_cum_dist(end:-1:1);
binom_probs = bnet_dist ./ bnet_reverse_cum_dist;  % Note that the last one is not needed ... 


binom_cumsum_vecs = zeros(length(bnet_dist)-1, chunk_size+1); % we do not need the last p=1

% Now generate the binom cumsum
for i=1:length(bnet_dist)-1
    binom_cumsum_vecs(i,:) =  binocdf([0:chunk_size],chunk_size,binom_probs(i));
end

% now remove all the not needed vectors which are outside the tol .. save
% for later .... 
