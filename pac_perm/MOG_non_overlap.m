% Test mixture-of-gaussians without overlsp
normed_labels = rand(1,nsamples) < p_labels;
normed_labels(1)=0; normed_labels(2)=1;  % Make sure we have both kinds
pos_ind = find(normed_labels); neg_ind = find(normed_labels==0);

normed_labels = normed_labels - mean(normed_labels);                
normed_labels = normed_labels ./ sqrt(normed_labels);


% Now randomize the data
normed_data = randn(nsamples, num_vars);

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
normed_labels = normed_labels ./ sqrt(normed_labels);


% Now randomize the data
normed_data = randn(nsamples, num_vars);

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
if(one_side_flag)
    [ val_corr ind_corr1] = sort(-corr_vec);
else
    [ val_corr ind_corr1] = sort(-abs(corr_vec));
end

top_vars = round(alpha * num_vars);

% Average over all the genes
inter_top = intersect(ind_corr1(1:top_vars),ind_corr(1:top_vars));
MOG_kept_frac_dist(i) = length(intersect(ind_corr1(1:top_vars),ind_corr(1:top_vars)));



