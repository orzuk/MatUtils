% Plot details of demographic model
%
% Input:
% D_cell - cell array with different demographic models
% index - representing the model chosen for each structure in D_cell
% log_like_mat - log-likelihood of each model
% k_vec - data (number of derived allele carriers)
% n_vec - data (number of individuals profiled)
% print_all_models  - flag (default 1)
%
function plot_time = demographic_model_plot(D_cell, index, log_like_mat, k_vec, n_vec, print_all_models)

ttt = cputime; 

AssignGeneralConstants; AssignRVASConstants;
plot_params.font_size=12;
if(~exist('print_all_models', 'var') || isempty(print_all_models))
    print_all_models =  1; % print all possible models
end
num_D = length(D_cell);
pseudo_count = 1;

figure; % subplot(1,2,1);  % Plot demographic models: pop. size vs. time
legend_vec = cell(num_D, 1);
for i=1:num_D
    N_vec = demographic_parameters_to_n_vec(D_cell{i}, index(i));
    semilogy(N_vec, 'linewidth', 2, 'color', color_vec(i)); hold on; % semilogy(N_vec_hat, 'r', 'linewidth', 2);
    legend_vec{i} = strdiff(D_cell{i}.name, 'Fitted.');
end
xlabel('Time (generations)'); ylabel('Population size');
if((num_D > 1) && print_all_models)
    title(['Best model: LL=' num2str(log_like_mat(index(2)), 5)]);
end
legend(legend_vec, 'fontsize', plot_params.font_size, 'location', 'southeast'); legend('boxoff');
%add_faint_grid(0.5);

if(~print_all_models)
    my_saveas(gcf, fullfile(exome_data_figs_dir, 'fitted_demographies'), {'epsc', 'pdf', 'jpg'}); % NEW: add .jpg for Robert
else
    good_inds = vec2row(find(abs(log_like_mat) < inf)) % filter to include only these 
    N_vec = demographic_parameters_to_n_vec(D_cell{1}, index(1)); % take first model
    sprintf('Computing models demographic distances ..')
    N_vec_cell = cell(D_cell{2}.num_params, 1); 
%    demographic_dist = zeros(length(good_inds), 1); 
    demographic_dist = zeros(D_cell{2}.num_params, 1); 
    demographic_dist_mat = zeros(length(good_inds)); % D_cell{2}.num_params);
    for i=1:D_cell{2}.num_params % good_inds % 1:D_cell{2}.num_params  % loop on all models of first
        N_vec_cell{i} = demographic_parameters_to_n_vec(D_cell{2}, i);
    end
    for i=1:D_cell{2}.num_params  % loop on all models of first
        demographic_dist(i) = demographic_models_distance(N_vec, N_vec_cell{i}); % N_vec_cell{good_inds(i)});  % find which one is most similar to the TRUE model
    end
    for i=1:length(good_inds) % 1:D_cell{2}.num_params  % loop on all models of first
        for j=(i+1):length(good_inds) % 2:D_cell{2}.num_params % compute all pairwise distances between models
            demographic_dist_mat(i,j) = demographic_models_distance(N_vec_cell{good_inds(i)}, N_vec_cell{good_inds(j)});   % find which ones are more similar to each other
        end
    end
    [~, min_I] = min(demographic_dist); % min_I = good_inds(min_I);  % find closest model with respect to demographic distance and print it
    semilogy(N_vec_cell{min_I}, 'linewidth', 2, 'color', color_vec(3), 'linestyle', '--'); %
    xlabel('Time (generations)'); ylabel('Population size');
    cur_title = get(gca, 'title');
    title([cur_title.String '; Similar model: LL=' num2str(log_like_mat(min_I), 5)]);
    legend([legend_vec' 'closest-demography'] , 'fontsize', plot_params.font_size); legend('boxoff');
    
    D_cell{end+1} = D_cell{end}; index(end+1) = min_I; D_cell{end}.index = min_I; % Add closest demography to true demography
    legend_vec{end+1} = 'closest-demography';
    
%    demographic_dist_mat = demographic_dist_mat+demographic_dist_mat'; % make symmetric
    demographic_embed = cmdscale(demographic_dist_mat+demographic_dist_mat'); % demographic_dist_mat(good_inds, good_inds)); % Run MDS
    add_faint_grid(0.5);
    my_saveas(gcf, fullfile(exome_data_figs_dir, 'fitted_demographies'), {'epsc', 'pdf', 'jpg'}); % NEW: add .jpg for Robert
    
    figure;  hold off; % PCA plot % subplot(1,2,2);
    scatter(demographic_embed(:,1), demographic_embed(:, 2), ...
        0.5*(log(max_cell(N_vec_cell(good_inds)))).^2.5, log_like_mat(good_inds)); hold on;
    plot(demographic_embed(good_inds == min_I,1), demographic_embed(find(good_inds == min_I), 2), 'kd'); colorbar; % closest model
    [~, best_I] = max(log_like_mat);
    plot(demographic_embed(find(good_inds == best_I),1), demographic_embed(find(good_inds == best_I), 2), 'r*'); colorbar; % highest-scoring model
    title('PCA LL of potential dem. models (size: max-N)');  xlabel('PC_1'); ylabel('PC_2');
    legend('All', 'Closest', 'Max-LL'); legend('boxoff');
end % print all models

plot_sfs = 1; % plot sfs for true and fitted model (the two should be close)
if(plot_sfs) % plot neutral sfs for demographic model
    rare_cumulative_per_gene = []; % set dummy variables
    target_size_by_class_vec = [mu_per_site, mu_per_site, mu_per_site]; % [neutral, null, missense]
    full_flag = 0; % use summary statistics
    null_w_vec = 1; % NULL_C % assume all alleles are 'null' (but s=0 so actually neutral)
    if(size(k_vec,2)==1) % duplicate data
        k_vec = repmat(k_vec, 1, length(D_cell));
    end
    if(size(n_vec,2)==1) % duplicate data
        n_vec = repmat(n_vec, 1, length(D_cell));
    end
    X = [vec2row(k_vec) vec2row(n_vec)]'; % Represent counts in a packed form
    n_sample = 200; % this should be determined in input
    LL_legend_vec = legend_vec;
    loglike_params = struct('null_w_vec', null_w_vec, 'include_phenotype', 0, ...
        'full_flag', full_flag, 'num_individuals', []);    % set parameters
    for i=1:length(D_cell)
        N_vec = demographic_parameters_to_n_vec(D_cell{i}, index(i));
        n_sample = min(n_sample, 2*N_vec(end-1));
        D_cell{i}.use_allele_counts=0;
        [log_like_mat_again{i}, P_poly_again{i}]  = ... % compute likelihood (here vary only alpha)
            compute_two_class_log_likelihood(0, 0, [], target_size_by_class_vec, D_cell{i}, ...
            X(:,i), [], [], loglike_params); % null_w_vec, 0, full_flag, []); % don't include phenotype !!
    end
    f_vec = k_vec ./ n_vec; % get empirical frequencies
    figure;  % Plot population SFS
    subplot(2,2,1);
    for i=1:length(D_cell)
        LL_legend_vec{i} = [legend_vec{i} ', LL=' num2str(round(log_like_mat_again{i}, 2))];
        N_vec = demographic_parameters_to_n_vec(D_cell{i}, index(i));
        if(~isfield(D_cell{i}, 'SFS'))
            [D_cell{i}.SFS.x_vec, D_cell{i}.SFS.p_vec, D_cell{i}.SFS.L, ~, k_vec_hat{i}, n_vec_hat{i}]  = ... % Compare neutral allele-freq distribution for different demographies
                compute_allele_freq_spectrum_from_demographic_model(D_cell{i}, 0, 'simulation', n_sample, mu_per_site); % simulate from neutral model
        else  % here compute only k_vec, n_vec
            max_num_alleles = 20000; num_alleles = max_num_alleles;
            p_vec_counts = round(D_cell{i}.SFS.p_vec(1,:) .* num_alleles);  allele_freq_vec = hist_to_vals(D_cell{i}.SFS.x_vec, p_vec_counts);
            k_vec_hat{i} = population_to_sample_allele_freq(allele_freq_vec, 2*N_vec(end-1), n_sample); % simulate a sample from population
            n_vec_hat{i} = repmat(n_sample, num_alleles, 1);
        end
        %        semilogx(x_vec ./ (2*N_vec(end-1)), p_vec); hold on;
        h_l(i) = semilogx(D_cell{i}.SFS.x_vec ./ (2*N_vec(end-1)), D_cell{i}.SFS.p_vec(1,:), color_vec(i), 'linewidth', 2); hold on;
        [unique_k_vec{i}, h_k_vec{i}] = unique_with_counts(k_vec(:,i));
    end
    for i=1:length(D_cell) % PLOT ALSO DATA
        semilogx(unique_k_vec{i} ./ max(n_vec(:,i)), h_k_vec{i} ./ length(k_vec(:,i)), [color_vec(i), '--'], 'linewidth', 2);
    end % loop on models
    xlabel('$f$ (population)', 'interpreter', 'latex'); ylabel('$\Psi(f)$ population', 'interpreter', 'latex');
    title('Population Frequencies', 'interpreter', 'latex');
    add_faint_grid(0.5, 0); set(gca, 'fontsize', plot_params.font_size-2);  x_tick = get(gca, 'XTick');
    set(gca, 'XTick', logspace(log10(x_tick(1)), log10(x_tick(end)), log10(x_tick(end)) - log10(x_tick(1)) +1)); % change ticks
    %    legend(h_l, LL_legend_vec, 'fontsize', plot_params.font_size); legend('boxoff');
    
    subplot(2,2,2);  % Plot sample SFS
    legend_with_LL = [LL_legend_vec'];
    for i=1:length(D_cell)
        [unique_k_vec_hat{i}, h_k_vec_hat{i}] = unique_with_counts(k_vec_hat{i});
        h_l(i) = semilogx(unique_k_vec_hat{i} ./ n_vec_hat{i}(1), h_k_vec_hat{i} ./ length(k_vec_hat{i}), ...
            color_vec(i), 'linewidth', 2); hold on;
        %        log_like_mat_direct(i) = sum(h_k_)vec_hat{i} .* log( h_k_vec_hat{i} ./ length(k_vec_hat{i}) ) );
        % New: Compute again likelihood for data given the ploted likelihoods:
        log_like_mat_direct(i) = relative_entropy_hist(unique_k_vec{i}, h_k_vec{i}, unique_k_vec_hat{i},  h_k_vec_hat{i}, 0, [1 1]); % compute LL. Add pseudo counts
        legend_with_LL{i} = [legend_with_LL{i} ' ; DLL=' num2str(log_like_mat_direct(i), 4)];
    end
    for i=1:length(D_cell) % plot data
        semilogx(unique_k_vec{i} ./ max(n_vec(:,i)), h_k_vec{i} ./ length(k_vec(:,i)), [color_vec(i), '--'], 'linewidth', 2);  % PLOT ALSO DATA
    end
    xlabel('$f$ (sample)', 'interpreter', 'latex'); ylabel('$\Psi(f)$ sample', 'interpreter', 'latex');
    legend(h_l, legend_with_LL, 'fontsize', plot_params.font_size-6, 'location', 'northwest'); legend('boxoff');
    title('Sample Frequencies', 'interpreter', 'latex');
    add_faint_grid(0.5, 0); set(gca, 'fontsize', plot_params.font_size-2);  x_tick = get(gca, 'XTick');
    set(gca, 'XTick', logspace(log10(x_tick(1)), log10(x_tick(end)), log10(x_tick(end)) - log10(x_tick(1)) +1)); % change ticks
    subplot(2,2,3);  hold off; % figure; % What is this figure?? is it empirical log-likelihood as function of derived allele count?? show also data on same plot!!!
    for i=1:length(D_cell) % loop on models. Here we consider polymorphic if smaller than MAX(n_vec) (should be modified to compare each allele to it's own n)
        h_l(i) = plot( unique(k_vec((k_vec(:,i)>0) & (k_vec(:,i) < max(n_vec(:,i))), i)), cumsum(P_poly_again{i}.LL_vec(2:end-1)), ['*' color_vec(i)]); hold on; % plot likelihood contribution for polymorphic alleles
    end
    for i=1:length(D_cell) % loop on models
        unique_k_vec_poly{i} = unique_k_vec{i}((unique_k_vec{i}>0) & (unique_k_vec{i}<max(n_vec(:,i))));
        h_k_vec_poly{i} = h_k_vec{i}((unique_k_vec{i}>0) & (unique_k_vec{i}<max(n_vec(:,i))));
        plot(unique_k_vec_poly{i}, cumsum((h_k_vec_poly{i}'+pseudo_count) .* log(P_poly_again{i}.sample_p_vec(2:end-1) + pseudo_count/length(k_vec(:,i)))), ['--' color_vec(i)]); hold on; % plot likelihood contribution for polymorphic alleles
    end
    %    plot(unique_k_vec, log(h_k_vec ./ length(k_vec)), 'm*', 'linewidth', 2);  % PLOT ALSO DATA
    xlabel('k (sample allele count.)'); ylabel('LogLike');
    add_faint_grid(0.5, 0); set(gca, 'fontsize', plot_params.font_size-2);
    %    legend(h_l, LL_legend_vec', 'fontsize', plot_params.font_size); legend('boxoff'); %  'Data']
    
    subplot(2,2,4); hold off; % plot counts ??
    legend_with_LL = [LL_legend_vec' 'Data'];
    for i=1:length(D_cell)
        h_l(i) = semilogx(unique_k_vec_poly{i} ./ max(n_vec(:,i)), P_poly_again{i}.sample_p_vec(2:end-1) ./ sum(P_poly_again{i}.sample_p_vec(2:end-1)), ...
            color_vec(i), 'linewidth', 2); hold on;
        log_like_mat_direct(i) = length(k_vec(:,i)) * relative_entropy_hist( ...
            unique_k_vec{i}, h_k_vec{i}, unique_k_vec{i},  P_poly_again{i}.sample_p_vec .* length(k_vec(:,i)), 0, [1 1]); % compute LL. Add pseudo counts
        legend_with_LL{i} = [legend_with_LL{i} ' ; DLL=' num2str(log_like_mat_direct(i), 4)]; % compute again
    end
    for i=1:length(D_cell)
        semilogx(unique_k_vec_poly{i} ./ max(n_vec(:,i)), h_k_vec_poly{i} ./ sum(h_k_vec_poly{i}), [color_vec(i), '--'], 'linewidth', 2);  % PLOT ALSO DATA  % / length(k_vec)
    end
    legend(h_l, legend_with_LL, 'fontsize', plot_params.font_size-6, 'location', 'northwest'); legend('boxoff');
    xlabel('k (sample allele count.)', 'interpreter', 'latex'); ylabel('$\psi(f)$ sample', 'interpreter', 'latex');
    add_faint_grid(0.5, 0); set(gca, 'fontsize', plot_params.font_size-2);  x_tick = get(gca, 'XTick');
    set(gca, 'XTick', logspace(log10(x_tick(1)), log10(x_tick(end)), log10(x_tick(end)) - log10(x_tick(1)) +1)); % change ticks
    my_saveas(gcf, fullfile(exome_data_figs_dir, 'fitted_sfs'), {'epsc', 'pdf', 'jpg'}); % NEW: add .jpg for Robert
    
    %     debug_SFS = 0;
    %     if(debug_SFS)
    %         x_vec1 = unique_k_vec_hat{i}(2:end-1) ./ n_vec_hat{i}(1);
    %         y_vec1 = h_k_vec_hat{i}(2:end-1) ./ length(k_vec_hat{i}); y_vec1 = y_vec1 ./ sum(y_vec1);
    %         x_vec2 = unique_k_vec_poly{i} ./ n_sample;
    %         y_vec2 = P_poly_again{i}.sample_p_vec(2:end-1) ./ sum(P_poly_again{i}.sample_p_vec(2:end-1));
    %         figure; semilogx(x_vec1, y_vec1); hold on;
    %         semilogx(x_vec2, y_vec2, 'r'); hold on;
    %         xlabel('f (allele freq.)'); ylabel('P(f) ???');
    %     end
end % if plot_sfs


% log_like_mat_again2 = ... % compute likelihood (here vary only alpha)
%     compute_two_class_log_likelihood(0, 0, [], target_size_by_class_vec, D_cell{2}, ...
%     X, [], [], null_w_vec, ...
%     0, full_flag, []); % don't include phenotype !!
%
% log_like_mat_again3 = ... % compute likelihood (here vary only alpha)
%     compute_two_class_log_likelihood(0, 0, [], target_size_by_class_vec, D_cell{3}, ...
%     X, [], [], null_w_vec, ...
%     0, full_flag, []); % don't include phenotype !!


% Now compare two models:
if(length(log_like_mat_again)>1)
    best_model_loglike = log_like_mat_again{2}
end
if(length(log_like_mat_again)>2)
    closest_demography_model_loglike = log_like_mat_again{3}
end

plot_time = cputime - ttt; 

