%  Generate one random sample from a Bayes net already given as a distribution (increased efficiency)
function sample_count = my_sample_bnet_generate_sample(bnet_cum_dist, num_samples, binom_probs,binom_cumsum_vecs, varargin)
% my_sample_bnet_generate_sample Generate a random sample from a Bayes net, (which is already given as a distribution) of size
% nsamples
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

% We now do it in the new way. Note that now the function returns the
% sufficient statistics, or the sample counts (NOT normalized by the number of samples taken), 
% and NOT the entire sample !!!!  This is for now
% good only for small BNETs. Thus we do not need to collect the samples later 



% First use the chunk_size as much times as we can 
sample_count = zeros(1, length(bnet_cum_dist));
chunk_size = size(binom_cumsum_vecs,2)-1;  % determine the chunk size
cur_num_samples = num_samples; % The initial number of samples 

for i=1:length(sample_count)-1 % go over the binomial distributions
    num_chunks =  floor(cur_num_samples/chunk_size);
    small_rand_vec = rand(1, num_chunks); 
    for j=1:num_chunks
        ind=bin_sear(binom_cumsum_vecs(i,:),small_rand_vec(j)); % Use binary search routine of Eran Ofek, to find the index of closest value (we need upper ..)        
        ind = max(ind,1); % avoid 0 index 

        if(binom_cumsum_vecs(i,ind) >= small_rand_vec(j))
            sample_count(i) = sample_count(i) + ind;
        else  % here our rand vec has passed the value so we have to add one 
            sample_count(i) = sample_count(i) + ind + 1;
% %             if   ( binom_cumsum_vecs(i,ind+1)<small_rand_vec(j))
% %                 errr_ind = ind+1
% %                 errrrr_binom_upper =     binom_cumsum_vecs(i,ind+1)
% %                 errr_rand_val = small_rand_vec(j)                
% %             end
        end
    end
    
    % Now we have to sample the remainder in the old way ... 
    if(   mod(cur_num_samples, chunk_size) > 0)
        rand_vec = rand(1,mod(cur_num_samples, chunk_size));

        sample_count(i) = sample_count(i) + sum (  rand_vec <= binom_probs(i)  ); 
%        mmmm = mod(cur_num_samples, chunk_size)
    end
    
    % Prepare for the next value
    cur_num_samples = cur_num_samples  - sample_count(i);
    
end

% Now update the last one 
sample_count(end) = num_samples-sum(sample_count);




% % % % Init to zero
% % % samples = zeros(1,num_samples)+length(bnet_cum_dist)-1; % We already start from the highest value!
% % % 
% % % % Now start sampling. For now we do it very simply : 
% % % % sample from this distribution
% % % 
% % % rand_vec = rand(1,num_samples);
% % % 
% % % % Do it by looping on values (not efficient) 
% % % for val=length(bnet_cum_dist)-1:-1:1  % We skip the first since its already updated  !!!!!!! 
% % %      samples(find(rand_vec <= bnet_cum_dist(val))) = val-1; 
% % % end



