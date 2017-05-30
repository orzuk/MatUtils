% Plot details of demographic model
%
% Input:
% D_cell - cell array with demographic models
% index - representing the correct model chosen
% log_like_mat - log-likelihood of each model
% k_vec - data (number of derived allele carriers)
% n_vec - data (number of individuals profiled)
%
function demographic_model_plot(D_cell, index, log_like_mat, k_vec, n_vec)

AssignGeneralConstants; AssignRVASConstants;
num_D = length(D_cell);


figure;  % Plot demographic models
legend_vec = cell(num_D, 1);
for i=1:num_D
    N_vec = demographic_parameters_to_n_vec(D_cell{i}, index(i));
    semilogy(N_vec, 'linewidth', 2, 'color', color_vec(i)); hold on; % semilogy(N_vec_hat, 'r', 'linewidth', 2);
    legend_vec{i} = D_cell{i}.name;
end
xlabel('Time (generations)'); ylabel('Population size'); title(['Best similar model: loglike=' num2str(log_like_mat(index(2)))]);
legend(legend_vec, 'fontsize', 14); legend('boxoff');

print_all_models = 1; % print all possible models
if(print_all_models)
    N_vec = demographic_parameters_to_n_vec(D_cell{1}, index(1));
    N_vec_cell = cell(D_cell{2}.num_params, 1); demographic_dist = zeros(D_cell{2}.num_params, 1);
    for i=1:D_cell{2}.num_params  % loop on all models of first
        N_vec_cell{i} = demographic_parameters_to_n_vec(D_cell{2}, i);
        
        % find which one is most similar to the true model
        demographic_dist(i) = demographic_models_distance(N_vec, N_vec_cell{i});
    end
    
    [min_demographic_dist, min_I] = min(demographic_dist);
    % Print the minimum with it's log-likelihood
    figure; semilogy(N_vec, 'linewidth', 2, 'color', color_vec(1)); hold on; % semilogy(N_vec_hat, 'r', 'linewidth', 2);
    semilogy(N_vec_cell{min_I}, 'linewidth', 2, 'color', color_vec(2));
    xlabel('Time (generations)'); ylabel('Population size'); title(['Most similar model: loglike=' num2str(log_like_mat(min_I))]);
    legend(legend_vec, 'fontsize', 14); legend('boxoff');
    
    D_cell{end+1} = D_cell{end}; index(end+1) = min_I; D_cell{end}.index = min_I; % Add closest demography to true demography
    legend_vec{end+1} = 'closest-demography';
end % print all models

plot_sfs = 1; % plot sfs for true and fitted model (Should be close)
if(plot_sfs) % plot neutral sfs for demographic model
    rare_cumulative_per_gene = []; % set dummy variables
    target_size_by_class_vec = [mu_per_site, mu_per_site, mu_per_site]; % [neutral, null, missense]
    full_flag = 0; % use summary statistics
    null_w_vec = 1; % NULL_C % assume all alleles are 'null' (but s=0 so actually neutral)
    X = [vec2row(k_vec) vec2row(n_vec)]'; % Represent counts in a packed form
    
    n_sample = 200; % this should be determined in input
    LL_legend_vec = legend_vec;
    
    for i=1:length(D_cell)
        D_cell{i}.use_allele_counts=0;
        [log_like_mat_again{i}, P_poly_again{i}]  = ... % compute likelihood (here vary only alpha)
            compute_two_class_log_likelihood(0, 0, [], rare_cumulative_per_gene, target_size_by_class_vec, D_cell{i}, ...
            X, [], [], null_w_vec, ...
            0, full_flag, []); % don't include phenotype !!
        LL_legend_vec{i} = [legend_vec{i} ', LL=' num2str(round(log_like_mat_again{i}, 2))];
    end
    
    figure;
    for i=1:length(D_cell)
        N_vec = demographic_parameters_to_n_vec(D_cell{i}, index(i));
        [x_vec_hat, p_vec_hat, L_correction_factor_hat, ~, k_vec_hat, n_vec_hat, weights_vec_hat]  = ... % Compare neutral allele-freq distribution for different demographies
            compute_allele_freq_spectrum_from_demographic_model(D_cell{i}, 0, 'simulation', n_sample, mu_per_site); % simulate from neutral model
        %        semilogx(x_vec ./ (2*N_vec(end-1)), p_vec); hold on;
        semilogx(x_vec_hat ./ (2*n_sample), p_vec_hat, color_vec(i)); hold on;
    end % loop on modelsai
    xlabel('f (allele freq.)'); ylabel('Psi(f)');
    legend(LL_legend_vec, 'fontsize', 14); legend('boxoff');
end % if plot_sfs


% log_like_mat_again2 = ... % compute likelihood (here vary only alpha)
%     compute_two_class_log_likelihood(0, 0, [], rare_cumulative_per_gene, target_size_by_class_vec, D_cell{2}, ...
%     X, [], [], null_w_vec, ...
%     0, full_flag, []); % don't include phenotype !!
%
% log_like_mat_again3 = ... % compute likelihood (here vary only alpha)
%     compute_two_class_log_likelihood(0, 0, [], rare_cumulative_per_gene, target_size_by_class_vec, D_cell{3}, ...
%     X, [], [], null_w_vec, ...
%     0, full_flag, []); % don't include phenotype !!


% Now compare:
best_model_loglike = log_like_mat_again{2}
closest_demography_model_loglike = log_like_mat_again{3}



% Compute how similar are two demographic models. Allow stretching of time (?)
% Input:
% N_vec1 - population size at each time for first model
% N_vec2 - population size at each time for second model
%
% Output:
% d - average difference in (log) population size between two models
%
function d = demographic_models_distance(N_vec1, N_vec2)

t1 = length(N_vec1); t2 = length(N_vec2);

if(t1>t2)
    d = demographic_models_distance(N_vec2, N_vec1);
else % here we know that t2>t1
    stretch_time=0;
    if(stretch_time) % compre all times stretched
        N_vec1 = interp1((1:t1).*t2./t1, N_vec1, (1:t2)', 'linear', N_vec1(1));
        d = mean((log(N_vec1)-log(N_vec2)).^2);
        
    else % compare only recent times
        d = mean((log(N_vec1)-log(N_vec2((t2-t1+1):end))).^2);
    end
end



