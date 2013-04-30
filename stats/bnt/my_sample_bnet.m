% Generate a random sample from a Bayes net - a modified implementation 
function samples = my_sample_bnet(bnet, num_samples, varargin)
% sample_bnet_nocell Generate a random sample from a Bayes net, of size
% nsamples
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

% set defauly params
nodes = length(bnet.dag);

% Init to zero
samples = zeros(1,num_samples);

%prrrr = 6

% Now start sampling. For now we do it very simply : 
% Transfer the BNT to probability distribution, and then sample from this
% distribution
bnet_dist = zuk_bnet_to_probs(bnet);

%tttttrrrr = 7


rand_vec = rand(1,num_samples);
bnet_cum_dist = cumsum(bnet_dist);


% Do it by looping on values (not efficient) 
for val=2^nodes:-1:1
     samples(find(rand_vec <= bnet_cum_dist(val))) = val-1; 
end



