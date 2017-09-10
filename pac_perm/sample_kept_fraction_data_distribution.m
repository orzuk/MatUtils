% Generate a sampling of the distribution for a given set of parameters
function [MOG_kept_frac_dist,kept_frac_dist] = sample_kept_fraction_data_distribution(rand_flag, one_side_flag,  true_corr_flag, 
true_corrs, nsamples, num_vars, alpha, num_iters, R, from_data_flag)
                                                              
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; FROM_DATA = 3; mix_GAUSSIAN=4;
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;

TRUE = 1; FALSE = 0;

kept_frac_dist = zeros(1, num_iters);

nsamples_inside_is = nsamples;
Nsamples=size(R.dat,2);
% Iterate sampling and determining the size of the kept fraction
genes_mu_vec = zeros(2,num_vars); genes_sigma_vec = zeros(2,num_vars);
all_pos_ind = find(R.Labels==1);
all_neg_ind = setdiff([1:Nsamples], all_pos_ind);

% Probability of a label to be one
p_labels = sum(R.Labels) / length(R.Labels);

genes_mu_vec(1,:) = mean(R.dat(:,all_pos_ind), 2);
genes_mu_vec(2,:) = mean(R.dat(:,all_neg_ind), 2);
genes_sigma_vec(1,:) = std(R.dat(:,all_pos_ind), [], 2);
genes_sigma_vec(2,:) = std(R.dat(:,all_neg_ind), [], 2);

% Here do the adjustment, to insure that we get the desired 'shrinked' correlations
mean_mu=tanh(true_corrs).*(genes_mu_vec(1,:)+genes_mu_vec(2,:))/2;
delta_mu=sqrt((p_labels*genes_sigma_vec(1,:).^2+...
    (1-p_labels)*genes_sigma_vec(2,:).^2)./(p_labels*(1-p_labels)*(1-(tanh(true_corrs)).^2)));
if(genes_mu_vec(1,:)>genes_mu_vec(2,:))
    genes_mu_vec(1,:)=mean_mu+delta_mu/2;
    genes_mu_vec(2,:)=mean_mu-delta_mu/2;
else
    genes_mu_vec(1,:)=mean_mu-delta_mu/2;
    genes_mu_vec(2,:)=mean_mu+delta_mu/2;
end
MOG_kept_frac_dist=[];
for i=1:num_iters
    if(mod(i,100) == 1)
        cur_iter = i
    end
    corr_vec=true_corrs;
    % Sample the true correlationsfrom the real distribution
    % Generate the sampling correlation vector
    sample_corr_vec = sample_noisy_correlations(true_corrs, nsamples, 1*GAUSSIAN); % This function should return the sample correlation

    % Now check how many of the top alpha were kept
    % count fraction of top genes kept
    if(one_side_flag)
        [ val_corr ind_corr] = sort(-sample_corr_vec);
    else
        [ val_corr ind_corr] = sort(-abs(sample_corr_vec));
    end
   
    
    just_before_if = 1;
    % here we have two samples !
    if(true_corr_flag == TWO_SAMPLED)

        here_inside_two_sampled = 2;
        % Here try to estimate the overlap from data
        %    figure; hist(sample_corr_vec, 300); title(['Corrleations from MODEL, ' num2str(nsamples)]); 
        corr_vec = sample_noisy_correlations(true_corrs, nsamples, 1*GAUSSIAN); % This function should return the sample correlation
        %         size([corr_vec ; sample_corr_vec])

    end

  
    if(one_side_flag)
         [ val_corr ind_corr1] = sort(-corr_vec);
    else
         [ val_corr ind_corr1] = sort(-abs(corr_vec));
    end
   top_vars = round(alpha * num_vars);

    kept_frac_dist(i) = length(intersect(ind_corr1(1:top_vars),ind_corr(1:top_vars)));
    MOG_non_overlap
end

