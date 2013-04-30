%  Generate many random samples from a Bayes net already given as a distribution (increased efficiency)
function sample_count = my_sample_bnet_generate_many_samples(bnet_cum_dist, num_samples, num_iters, binom_probs,binom_cumsum_vecs, varargin)
% zuk_sample_bnet_generate_sample Generate a random sample from a Bayes net, (which is already given as a distribution) of size
% nsamples
%
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


% First use the chunk_size as much times as we can
sample_count = zeros(num_iters, length(bnet_cum_dist));
chunk_size = size(binom_cumsum_vecs,2)-1;  % determine the chunk size


for iter = 1:num_iters
    cur_num_samples = num_samples; % The initial number of samples

    for i=1:length(bnet_cum_dist)-1 % go over the binomial distributions
        num_chunks =  floor(cur_num_samples/chunk_size);
        small_rand_vec = rand(1, num_chunks);
        for j=1:num_chunks
            % Use binary search routine of Eran Ofek, to find the index of
            % closest value (we need upper ..)
            ind=bin_sear(binom_cumsum_vecs(i,:),small_rand_vec(j));
            ind = max(ind,1); % avoid 0 index

            if(binom_cumsum_vecs(i,ind) >= small_rand_vec(j))
                sample_count(iter,i) = sample_count(iter,i) + ind;
            else  % here our rand vec has passed the value so we have to add one
                sample_count(iter,i) = sample_count(iter,i) + ind + 1;
            end
        end

        % Now we have to sample the remainder in the old way ...
        if(   mod(cur_num_samples, chunk_size) > 0)
            
%             cur_num_samples
%             chunk_size
%             do_some_leftovers = 999999999
            
            rand_vec = rand(1,mod(cur_num_samples, chunk_size));

            sample_count(iter,i) = sample_count(iter,i) + sum (  rand_vec <= binom_probs(i)  );
        end

        % Prepare for the next value
        cur_num_samples = cur_num_samples  - sample_count(iter,i);

    end

    % Now update the last one
    sample_count(iter,end) = num_samples-sum(sample_count(iter,:));

end



