% Generate a sampling of the distribution for a given set of parameters
% New: We added as a 4th argument, the 'true' Z values. 
function [ kept_frac_dist f_all_genes f_all_genes_one_NTOP] = sample_kept_fraction_distribution(rand_flag, one_side_flag,  true_corr_flag, true_Zs, ...
    dist_std, nsamples, num_vars, alpha, num_iters, R, from_data_flag,dist_mean,prior, sig_sig, coeffs_vec, corr_mat, sampling_flag)

global BOOTSTRAP NON_OVERLAP MOG_NON_OVERLAP;
global UNIFORM GAUSSIAN LINEAR FROM_DATA mix_GAUSSIAN student_t;
global ONE_SIDE TWO_SIDES;
global TRUE_AND_SAMPLED TWO_SAMPLED;
global TRUE FALSE;


Ngenes = size(R.dat, 1)
Nsamples = size(R.dat, 2)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       NEW: For each gene fit its Gaussian parameters for good and poor rognosis
genes_mu_vec = zeros(2,Ngenes); genes_sigma_vec = zeros(2,Ngenes);
all_pos_ind = find(R.Labels==1);
all_neg_ind = setdiff([1:Nsamples], all_pos_ind);

% Probability of a label to be one
p_labels = sum(R.Labels) / length(R.Labels);

genes_mu_vec(1,:) = mean(R.dat(:,all_pos_ind), 2);
genes_mu_vec(2,:) = mean(R.dat(:,all_neg_ind), 2);
genes_sigma_vec(1,:) = std(R.dat(:,all_pos_ind), [], 2);
genes_sigma_vec(2,:) = std(R.dat(:,all_neg_ind), [], 2);

% Do the correction to get appropriate correlations
% Here do the adjustment, to insure that we get the desired 'shrinked' correlations
mean_mu=(genes_mu_vec(1,:)+genes_mu_vec(2,:))/2;
delta_mu=tanh(true_Zs).*sqrt((p_labels*genes_sigma_vec(1,:).^2+...
    (1-p_labels)*genes_sigma_vec(2,:).^2)./(p_labels*(1-p_labels)*(1-(tanh(true_Zs)).^2)));
if(genes_mu_vec(1,:)>genes_mu_vec(2,:))
    genes_mu_vec(1,:)=mean_mu+delta_mu/2;
    genes_mu_vec(2,:)=mean_mu-delta_mu/2;
else
    genes_mu_vec(1,:)=mean_mu-delta_mu/2;
    genes_mu_vec(2,:)=mean_mu+delta_mu/2;
end
% % % % % % figure; hold on; hist(genes_mu_vec(1,:), 100); title('pos mu labels hist');
% % % % % % figure; hold on; hist(genes_mu_vec(2,:), 100); title('neg mu labels hist');
% % % % % % figure; hold on; hist(genes_sigma_vec(1,:), 100); title('pos sigma labels hist');
% % % % % % figure; hold on; hist(genes_sigma_vec(2,:), 100); title('neg sigma labels hist');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



kept_frac_dist = zeros(1, num_iters);

nsamples_inside_is = nsamples;
num_of_Gaussians=length(dist_std);


% Take the sqrt of a matrix
corr_mat_sqrt = [];


% We now compute the overlap of each gene seperately to see what will
% happen
f_all_genes = zeros(2,num_vars);  % We have 00 01 10 11 , i.e. we have redundancy
f_all_genes_one_NTOP = zeros(2, num_vars);

i=1;

% Probability of a label to be one
p_labels = sum(R.Labels) / length(R.Labels);


% Iterate sampling and determining the size of the kept fraction
while(i <= num_iters)

    update_i=1;
    if(mod(i,100) == 1)
        cur_iter = i
    end

    if(from_data_flag == FALSE)
        % Sample the true correlations
        corr_vec = sample_true_correlations(rand_flag, num_vars, dist_mean, dist_std,prior,num_of_Gaussians); % Take zero mean and one variance
        % Generate the sampling correlation vector
        sample_corr_vec = sample_noisy_correlations(corr_vec, nsamples, 1*GAUSSIAN, sig_sig, coeffs_vec, corr_mat_sqrt); % This function should return the sample correlation

        % Now check how many of the top alpha were kept
        % count fraction of top genes kept
        if(one_side_flag)
            [ val_corr ind_corr] = sort(-sample_corr_vec);
        else
            [ val_corr ind_corr] = sort(-abs(sample_corr_vec));
        end
    end

    % here we have two samples !
    if(true_corr_flag == TWO_SAMPLED)

        % Here try to estimate the overlap from data
        if(from_data_flag == TRUE)
            
            % New: allow for MoG sampling from the data
            if(sampling_flag == MOG_NON_OVERLAP)
                
                % first randomise the labels
                normed_labels = rand(1,nsamples) < p_labels;
                normed_labels(1)=0; normed_labels(2)=1;  % Make sure we have both kinds
                pos_ind = find(normed_labels); neg_ind = find(normed_labels==0);
                
                normed_labels = normed_labels - mean(normed_labels);                
                normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
                
                
                % Now randomize the data
                normed_data = randn(nsamples, Ngenes);
               
                normed_data(pos_ind,:) = normed_data(pos_ind,:) .* repmat(genes_sigma_vec(1,:), length(pos_ind),1);
                normed_data(pos_ind,:) = normed_data(pos_ind,:) +  repmat(genes_mu_vec(1,:),length(pos_ind),1);
                normed_data(neg_ind,:) = normed_data(neg_ind,:) .* repmat(genes_sigma_vec(2,:),length(neg_ind),1);
                normed_data(neg_ind,:) = normed_data(neg_ind,:) +  repmat(genes_mu_vec(2,:),length(neg_ind),1);
           
                % Finally normalize the data
                normed_data = normed_data - repmat(mean(normed_data, 1), nsamples, 1);
                normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 1)), nsamples, 1);
                
                sample_corr_vec = normed_data' * normed_labels';
                
                
                % Now do everything again
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                 % first randomise the labels
                normed_labels = rand(1,nsamples) > p_labels;
                normed_labels(1)=0; normed_labels(2)=1;  % Make sure we have both kinds
                pos_ind = find(normed_labels); neg_ind = find(normed_labels==0);
                
                normed_labels = normed_labels - mean(normed_labels);                
                normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
                
                
                % Now randomize the data
                normed_data = randn(nsamples, Ngenes);
               
                normed_data(pos_ind,:) = normed_data(pos_ind,:) .* repmat(genes_sigma_vec(1,:), length(pos_ind),1);
                normed_data(pos_ind,:) = normed_data(pos_ind,:) +  repmat(genes_mu_vec(1,:),length(pos_ind),1);
                normed_data(neg_ind,:) = normed_data(neg_ind,:) .* repmat(genes_sigma_vec(2,:),length(neg_ind),1);
                normed_data(neg_ind,:) = normed_data(neg_ind,:) +  repmat(genes_mu_vec(2,:),length(neg_ind),1);
           
                % Finally normalize the data
                normed_data = normed_data - repmat(mean(normed_data, 1), nsamples, 1);
                normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 1)), nsamples, 1);
                
                corr_vec = normed_data' * normed_labels';
                
                if(one_side_flag)
                    [ val_corr ind_corr] = sort(-sample_corr_vec);
                else
                    [ val_corr ind_corr] = sort(-abs(sample_corr_vec));
                end

            else % Here do the old way - really from the data
                cur_perm = randperm(length(R.Labels));

                % Calculate the mean correlation of each gene with survival
                normed_labels = R.Labels(cur_perm([1:nsamples])) - mean(R.Labels(cur_perm([1:nsamples])));

                % We need to avoid a case where all the flags or the
                % data's!!! are equal, so we can't define a correlation
                if(min(normed_labels) == max(normed_labels))
                    update_i = 0;
                end

                normed_data = R.dat(:,cur_perm([1:nsamples])); %  - repmat(mean(R.dat(:,cur_perm([1:nsamples])), 2), 1, nsamples);
                normed_data = normed_data - repmat(mean(normed_data, 2), 1, nsamples);

                if(~isempty(find(min(normed_data,[],2) == max(normed_data,[],2)))) % We simultaniously check for all rows
                    update_i = 0;
                end

                if(update_i == 1)
                    normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
                    normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, nsamples);
                    %sample_corr_vec = sum(normed_data .* repmat(normed_labels, num_vars, 1),2);
                    sample_corr_vec = normed_data * normed_labels';

                    % Calculate the mean correlation of each gene with survival
                    normed_labels = R.Labels(cur_perm([1*nsamples+1:2*nsamples])) - mean(R.Labels(cur_perm([1*nsamples+1:2*nsamples])));
                    normed_labels = normed_labels ./ sqrt(sum(normed_labels.^2));
                    normed_data = R.dat(:,cur_perm([1*nsamples+1:2*nsamples])); % - repmat(mean(R.dat(:,cur_perm([1*nsamples+1:2*nsamples])), 2), 1, nsamples);
                    normed_data = normed_data - repmat(mean(normed_data, 2), 1, nsamples);
                    normed_data = normed_data ./ repmat(sqrt(sum(normed_data.^2, 2)),1, nsamples);
                    %corr_vec = sum(normed_data .* repmat(normed_labels, num_vars, 1),2);
                    corr_vec = normed_data * normed_labels';

                    if(one_side_flag)
                        [ val_corr ind_corr] = sort(-sample_corr_vec);
                    else
                        [ val_corr ind_corr] = sort(-abs(sample_corr_vec));
                    end
                end
            end % if MOG_NON_OVERLAP

        else % TRUE   Here data-flag is false. Sample from model again !
            

            save_corr_vec = corr_vec; % Save the previous correlations
            corr_vec = sample_noisy_correlations(corr_vec, nsamples, 1*GAUSSIAN, sig_sig, coeffs_vec, corr_mat_sqrt); % This function should return the sample correlation


            %         size([corr_vec ; sample_corr_vec])
            mean_corrs_from_model = mean(sample_corr_vec);
            std_corrs_from_model = std(sample_corr_vec);
            nada_ind_corr_here = 0;
        end

    end
%    sort_and_inter_ttt = cputime;
    

    if(update_i == 1)
        if(one_side_flag)
            [ val_corr ind_corr1] = sort(-corr_vec);
        else
            [ val_corr ind_corr1] = sort(-abs(corr_vec));
        end
        
        top_vars = round(alpha * num_vars);

        % Average over all the genes
        inter_top = intersect(ind_corr1(1:top_vars),ind_corr(1:top_vars));
        inter_bottom = intersect(ind_corr1(top_vars+1:end),ind_corr(top_vars+1:end));
        kept_frac_dist(i) = length(intersect(ind_corr1(1:top_vars),ind_corr(1:top_vars)));
        
        
        % Here deal with each gene seperately
        f_all_genes(1,inter_top) = f_all_genes(1,inter_top)+1;
        f_all_genes(2,inter_bottom) = f_all_genes(1,inter_bottom)+1;
        
        f_all_genes_one_NTOP(1,ind_corr(1:top_vars)) = f_all_genes_one_NTOP(1,ind_corr(1:top_vars))+1;
        f_all_genes_one_NTOP(2,ind_corr(top_vars+1:end)) = f_all_genes_one_NTOP(2,ind_corr(top_vars+1:end))+1;
        
    end
%    sort_and_inter_ttt = cputime - sort_and_inter_ttt 
    
    i=i+update_i;

end

f_all_genes = f_all_genes./ num_iters;
f_all_genes_one_NTOP = f_all_genes_one_NTOP ./ num_iters; f_all_genes_one_NTOP = f_all_genes_one_NTOP .^ 2;

% Output only the prob. of being in N_top in one group 
% given that I'm in N_top in other group

%%%%% Currently we do not calculate the ratio:
%%%%% f_all_genes = 2 .* f_all_genes(1,:) ./ (1 + f_all_genes(1,:) + f_all_genes(2,:)); 

% if(from_data_flag == FALSE)
%     MEAN_IS = mean(save_corr_vec)
%     STD_IS = std(save_corr_vec)
%     figure; hist(save_corr_vec, 100); title('Sampled the correlations'); xlabel(['Mean is ' num2str(MEAN_IS) ' std is ' num2str(STD_IS)]);
% end



