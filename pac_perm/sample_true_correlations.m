% Now function sample_true_correlations
% The function samples a vector of original (true) corrleations, according to a given distribution
function samp_corr_vec = sample_true_correlations(rand_flag, num_vars, dist_mean, dist_std,prior,num_of_Gaussians)

UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; mix_GAUSSIAN=4;

% Sample the true correlations
if(rand_flag == UNIFORM)
    samp_corr_vec = (rand(1, num_vars)-0.5)*2*sqrt(3)*dist_std + dist_mean; 
end
if(rand_flag == GAUSSIAN)
    % gaussian
    samp_corr_vec = randn(1, num_vars)*dist_std+dist_mean;   
end
if(rand_flag == mix_GAUSSIAN)
    % gaussian
    samp_corr_vec=[];
    for i=1:num_of_Gaussians-1
        n(i)=round(prior(i)*num_vars);
        samp_corr_vec(end+1:end+n(i)) = randn(1, n(i))*dist_std(i)+dist_mean(i);   
    end
    n(end+1)=num_vars-sum(n);
    samp_corr_vec(end+1:end+n(end)) = randn(1, n(end))*dist_std(end)+dist_mean(end);   
end

if(rand_flag == LINEAR)
    % linear triangle-shaped correlations distribution
    helper_vec = rand(2,num_vars);  
    flip_index_vec = helper_vec(1,:) > helper_vec(2,:);
    helper_vec(1,flip_index_vec) = helper_vec(1,flip_index_vec)-1;
    
    samp_corr_vec = helper_vec(1,:)*sqrt(6)*dist_std+dist_mean  ;   
end


samp_corr_vec = sort(samp_corr_vec);

