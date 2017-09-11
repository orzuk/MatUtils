% Simulate parents and offspring and regress them for the LP model
%
% Input:
% iters - number of individuals to simulate
% N - number of pathways
% h_x - heritability of each pathway
% h_shared_env - shared environment componenet
% plot_flag - plot figures or not
% s_vec - vector of selection threshols 
% 
% Output:
% beta - regression coefficient of response on selection strength
% h_all - heritability explained by all loci (bottom-up)
% midparent_vec - vector of possible deviation for the mid-parent value
% offspring_mean_vec - vector of mean offspring deviation for each mid-parent deviation
% mean_selected_parents_phenotype - mean mid-parents phenotypes after selection
% mean_selected_offspring_phenotype - mean offspring phenotypes after selection 
%
function [beta h_all midparent_vec offspring_mean_vec ...
    mean_selected_parents_phenotype mean_selected_offspring_phenotype ...
    mid_parent_qtl_vec offspring_qtl_vec selected_inds ...
    mu_s_vec mu mu_offspring_vec mu_offspring] = ...
    simulate_directional_selection(iters, N, h_x, h_shared_env, ...
    plot_flag, s_vec)


figs_dir = '../../common_disease_model/docs/pnas/genetic_interactions/figs/breeding';
if(~exist('h_shared_env', 'var') || isempty(h_shared_env))
    h_shared_env = 0;
end
if(~exist('plot_flag', 'var') || isempty(plot_flag))
    plot_flag = 1;
end

compute_type='sampling'; % 'numeric'; % 'analytic';
selection_mode='lowerbound';

[mu_s sigma_s] = maxnormstat(N); 
h_unique_env = 1-h_x-h_shared_env;
parents_shared_env_vec = randn(N, iters) .* sqrt(h_shared_env);

paternal_genetic_vec = randn(2, N, iters) .* sqrt(h_x) ./ sqrt(2); % split genotype of parents to two parts
paternal_unique_env_vec = randn(N, iters) .* sqrt(h_unique_env);
paternal_env_vec = parents_shared_env_vec + paternal_unique_env_vec;
paternal_qtl_vec = ...
    max(paternal_env_vec + reshape(sum(paternal_genetic_vec), N, iters),[],1);
maternal_genetic_vec = randn(2, N, iters) .* sqrt(h_x) ./ sqrt(2);
maternal_unique_env_vec = randn(N, iters) .* sqrt(h_unique_env);
maternal_env_vec = parents_shared_env_vec + maternal_unique_env_vec;
maternal_qtl_vec = ...
    max(maternal_env_vec + reshape(sum(maternal_genetic_vec), N, iters),[],1);

maternal_qtl_vec = (maternal_qtl_vec - mu_s) ./ sigma_s;
paternal_qtl_vec = (paternal_qtl_vec - mu_s) ./ sigma_s;

mid_parent_qtl_vec = 0.5*(paternal_qtl_vec + maternal_qtl_vec);

paternal_choose_vec = rand(N, iters) > 0.5; % choose which one to pass on from paternal genome
maternal_choose_vec = rand(N, iters) > 0.5; % choose which one to pass on from maternal genome
offspring_genetic_vec = zeros(2, N, iters);
offspring_genetic_vec(1,:,:) = ...
    paternal_genetic_vec(1,:,:) .* reshape(paternal_choose_vec, 1, N, iters) + ...
    paternal_genetic_vec(2,:,:) .* (1-reshape(paternal_choose_vec, 1, N, iters));
offspring_genetic_vec(2,:,:) = ...
    maternal_genetic_vec(1,:,:) .* reshape(maternal_choose_vec, 1, N, iters) + ...
    maternal_genetic_vec(2,:,:) .* (1-reshape(maternal_choose_vec, 1, N, iters));

%offspring_shared_env_vec = sum(offspring_shared_env_vec); 
% maternal_shared_env_vec_chosen = maternal_shared_env_vec(1,:,:) .* 
% offspring_shared_env_vec = (maternal_shared_env_vec + paternal_shared_env_vec) ./ 2; % Take half from each parent % sqrt(2); 
%h_offspring_unique_env = 1-h_x-h_shared_env;

offspring_unique_env_vec = randn(N, iters) .* sqrt(h_unique_env); % This shouls be split to shared and unique environment
offspring_env_vec = parents_shared_env_vec + offspring_unique_env_vec;

offspring_qtl_vec = ...  % generate phenotype for offspring of two individuals
    max(offspring_env_vec + reshape(sum(offspring_genetic_vec), N, iters),[],1);


mz_genetic_vec = paternal_genetic_vec(1,:,:) + paternal_genetic_vec(2,:,:);
mz_qtl_vec = ... % generate also a phenotype for and mz twin of an individual
    max(offspring_env_vec +  reshape(mz_genetic_vec, N, iters),[],1);
offspring_qtl_vec = (offspring_qtl_vec - mu_s) ./ sigma_s;
mz_qtl_vec = (mz_qtl_vec - mu_s) ./ sigma_s;


mu = mean(mid_parent_qtl_vec); % mean of population
mu_offspring = mean(offspring_qtl_vec);
mu_mz = mean(mz_qtl_vec)
if(~exist('s_vec', 'var') || isempty(s_vec))
    s_vec = [-5:0.05:1.2]; % allow only people above this value to mate
%    s_vec = [-5:2:1.1]; % allow only people above this value to mate
end
s_vec_inds=ones(length(s_vec),1);
mu_s_vec = zeros(length(s_vec), 1); % differences in phenotypic value for parents
mu_offspring_vec  = zeros(length(s_vec), 1); % differences in phenotypic value for offspring
mu_offspring_vec_numeric  = zeros(length(s_vec), 1); % differences in phenotypic value for offspring
mu_mz_vec  = zeros(length(s_vec), 1); % differences in phenotypic value for offspring


for i=1:length(s_vec) % loop on different thresholds for selection
    sprintf('run selection threshold %ld of %ld', i, length(s_vec))
    switch selection_mode
        case {'right', 'lowerbound'} % allow only individuals with phenotype above this value to reproduce
            good_inds = find(min(paternal_qtl_vec, maternal_qtl_vec) > s_vec(i));
        case 'point' % allow only individuals with exact phenotype to reproduce
            good_inds = intersect(find(min(paternal_qtl_vec, maternal_qtl_vec) > s_vec(i)), ...
                find(max(paternal_qtl_vec, maternal_qtl_vec) < s_vec(i)+0.01));
            if(isempty(good_inds))
                s_vec_inds(i)=0;
            end
    end
    switch compute_type
        case 'sampling'
            mu_s_vec(i) = mean(mid_parent_qtl_vec(good_inds));
        case 'numeric'
            mu_s_vec(i) =  quadl(@(x)x.*maxnormpdf(x,N), s_vec(i), mu_s+5*sigma_s) ./ ...
                (1-maxnormcdf(s_vec(i),N)); % this is conditional mean assuming we're above s_vec(i)
            mu_offspring_vec_numeric(i) = quadl(@(s) 1 - ...  % compute offspring shift as if it's MZ twin
                conditional_offspring_cumulative_LP(s, h_x, s_vec(i)).^N - ...
                conditional_offspring_cumulative_LP(-s, h_x, s_vec(i)).^N, 0, 10*sigma_s);
            %            mu_offspring_vec_numeric(i) = (mu - normcdf(s_vec(i)).^N * mu_offspring_vec_numeric(i)) / ...
            %                (1-normcdf(s_vec(i)))^N; % take the opposite selection direction
    end
    mu_offspring_vec(i) = mean(offspring_qtl_vec(good_inds));  % this we know how to do only numerically ..
    mu_mz_vec(i) = mean(mz_qtl_vec(good_inds));  % this we know how to do only numerically ..
end
mu_is = mu, mu_analytic_is = mu_s
[mu_again sigma k_of_N_hist k_of_N_bins_loc corr_vec ...
    h_all h_pop h_x_vec qtl_R same_inds_prob] = ...
    compute_k_of_N_gaussian_statistics([0 0], [1 1], h_x, h_shared_env, [], ...
    'MAX', [], ...
    N, 10, 'numeric', 0, {'ACE', 'ADE', 'PO'});
h_all_is = h_all
h_pop



beta = polyfit(mid_parent_qtl_vec, offspring_qtl_vec, 1);
beta = beta(1)


x_vec = min(mu_s_vec-mu-0.1):0.01:max(mu_s_vec-mu+0.1); % 0:0.1:1;
beta_offspring = polyfit(mu_s_vec(find(s_vec_inds)) - mu, mu_offspring_vec(find(s_vec_inds)) - mu_offspring, 1);

if(plot_flag) % see if we want to plot
    figure; hold on; % Plot regression of delta_mu and delta_s for different selection coefficients.
    plot(mu_s_vec - mu, mu_offspring_vec - mu_offspring, '.', 'markersize', 10);
    %beta_offspring_normalized = corr(vec2column(mu_s_vec - mu), ...
    %    vec2column(mu_offspring_vec - mu_offspring));
    xlabel('\Delta s', 'fontsize', 14, 'fontweight', 'bold'); %  = \mu_s - \mu');
    ylabel('\Delta \mu', 'fontsize', 14, 'fontweight', 'bold'); %  = \mu_{offspring} - \mu');
    title([repmat(' ', 1, 80) 'supp-fig8'], 'fontsize', 16, 'fontweight', 'bold');
    % % % % % %   No title!!!
    % % % % %     title(['Selection vs. difference in offspring value. LP(' ...
    % % % % %         num2str(N) ', ' num2str(h_x*100,2) '%). \beta=' num2str(beta_offspring(1)*100,3) '%. h_{loci}='  ...
    % % % % %         num2str(h_all*100,3) '%  h_{pop}^2(PO)=' ...
    % % % % %         num2str(h_pop{3}*100,3) '%']);
    
    
    % % %         '%  h_{pop,DZ}^2=' num2str(2*qtl_R(2)*100,3) ...
    % % %         '% h_{pop,MZ}=' num2str(qtl_R(1)*100,3) '%']);
    plot(x_vec, x_vec .* beta_offspring(1) + beta_offspring(2), 'r', 'linewidth', 2);
    xlim([0, max(mu_s_vec-mu).*1.01]);
    ylim([0, max(max(mu_offspring_vec - mu_offspring).*1.01, ...
        max(max(mu_s_vec-mu).*1.01 .* beta_offspring(1) + beta_offspring(2)).*1.01)])
    switch compute_type
        case 'numeric'
            plot(mu_s_vec - mu, mu_offspring_vec_numeric - mu_offspring, 'g', 'markersize', 10);
    end
%    plot(mu_s_vec - mu, mu_mz_vec - mu_mz, '.m', 'markersize', 10); % plot empirical mz difference
%    legend({'offspring (samp.)', 'linear-fit', '\Delta mz (numeric.)', '\Delta mz (samp.)'}, 2);
%    legend('boxoff');
    my_saveas(gcf, fullfile(figs_dir, ...
        'response_to_selection'), {'epsc', 'pdf'});
end % if plot
%    beta_offspring = beta_offspring(1);


paternal_qtl_vec_normalized = NormalizeData(paternal_qtl_vec);
maternal_qtl_vec_normalized = NormalizeData(maternal_qtl_vec);
mid_parent_qtl_vec_normalized = NormalizeData(mid_parent_qtl_vec);
offspring_qtl_vec_normalized = NormalizeData(offspring_qtl_vec);
rho = corr(mid_parent_qtl_vec', offspring_qtl_vec')


num_bins = 1000; % increase number of bins 
[mid_parent_hist x_vec] = hist(mid_parent_qtl_vec, num_bins); % perform hist
bin_size = x_vec(2)-x_vec(1);
mid_parent_qtl_vec_binned = round(mid_parent_qtl_vec ./ bin_size); % .*bin_size;
mid_parent_qtl_vec_binned = mid_parent_qtl_vec_binned-min(mid_parent_qtl_vec_binned)+1;
offspring_hist = ...
    accumarray(mid_parent_qtl_vec_binned', offspring_qtl_vec', [], @mean);
offspring_counts = ...
    accumarray(mid_parent_qtl_vec_binned', ones(iters,1), [], @sum);
good_bins = find(offspring_counts > 100); % throw away sparse bins

rand_inds = randperm(iters); rand_inds = rand_inds(1:min(iters, 10000)); % don't plot everything ...
[mid_parent_qtl_vec_normalized_sorted sort_perm] = sort(mid_parent_qtl_vec_normalized);
offspring_qtl_vec_normalized_sorted = smooth(offspring_qtl_vec_normalized(sort_perm), 10);
paternal_qtl_vec_normalized_sorted = paternal_qtl_vec_normalized(sort_perm);
maternal_qtl_vec_normalized_sorted = maternal_qtl_vec_normalized(sort_perm);

offspring_qtl_vec_sorted = offspring_qtl_vec(sort_perm);
mid_parent_qtl_vec_sorted = mid_parent_qtl_vec(sort_perm); 


T = 1; % selected_inds = find(mid_parent_qtl_vec_normalized > T);
selected_inds = find(min(paternal_qtl_vec_normalized, maternal_qtl_vec_normalized) > T);
selected_rand_inds = intersect(rand_inds, selected_inds); % choose only the ones above threshold 
mean_selected_parents_phenotype = mean(mid_parent_qtl_vec(selected_inds)); % no need to normalize again!
mean_selected_offspring_phenotype = mean(offspring_qtl_vec(selected_inds)); % no need to normalize again! 
h_inbreeding = sqrt(2) * mean_selected_offspring_phenotype / mean_selected_parents_phenotype; % Delta mu over delta s

if(plot_flag)
    figure; hold on; % Plot regression of delta_mu and delta_s for different selection coefficients.
    plot(mid_parent_qtl_vec(rand_inds), offspring_qtl_vec(rand_inds), '.'); % all pointts
    plot(mid_parent_qtl_vec(selected_rand_inds), offspring_qtl_vec(selected_rand_inds), 'r.'); % selected pointts
    plot(x_vec(good_bins), offspring_hist(good_bins), 'g', 'linewidth', 2.5);     % moving average 
    plot(mid_parent_qtl_vec_sorted, smooth(offspring_qtl_vec_sorted, 5000), 'c')
    xlabel('mid-parent phenotype'); ylabel('offspring phenotype');
    beta_fit = polyfit(mid_parent_qtl_vec, offspring_qtl_vec, 1);
    plot(x_vec, x_vec .*beta_fit(1)+beta_fit(2), 'r', 'linewidth', 2);
    plot(x_vec, x_vec .* h_all, '--m', 'linewidth', 2); % plot line with h_all
    line([T T], [-2 2], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
    plot(mean_selected_parents_phenotype, mean_selected_offspring_phenotype, 'xk', 'linewidth', 10);
    %    text(mean_selected_parents_phenotype, mean_selected_offspring_phenotype, ...
    %        'x', 'fontsize', 20, 'fontweight', 'bold');
    legend({'data', 'selected', 'local-mean', 'smoothed', 'linear-regression', ...
        'y=h_{all}^2 x', 'selection-threshold'}, 2); legend('boxoff');
    model_str = ['LP(' num2str(N) ', ' num2str(h_x*100) '%)'];
    %    title(['Regression \beta=' num2str(beta_fit(1),3) ', h_{pop}^2(PO)=' num2str(h_pop{3}*100,3) ...
    %        '% r_{MZ}=' num2str(qtl_R(1),3) ' r_{DZ}=' num2str(qtl_R(2),3) ...
    %        ', h_{all}^2=' num2str(h_all*100,3) '%']);
    title([model_str ', h_{pop}^2(PO)=' num2str(h_pop{3}*100,3) ...
        '%, h_{\Delta \mu/\Delta S}^2=' num2str(h_inbreeding*100,3) ...
        '%, h_{all}^2=' num2str(h_all*100,3) '%']);
%    ylim([-2 2]);
    my_saveas(gcf, fullfile(figs_dir, 'response_to_selection'), 'epsc');
    
    
    
    x_mat = repmat(x_vec, length(x_vec), 1);
    offspring_hist = offspring_hist(1:end-1);
    offspring_hist_mat = repmat(offspring_hist, 1, length(offspring_hist));
    figure; hold on; % Plot everything in 3-D
    plot3(paternal_qtl_vec_normalized_sorted(rand_inds), maternal_qtl_vec_normalized_sorted(rand_inds), ...
        offspring_qtl_vec_normalized_sorted(rand_inds), '.');
    surf(x_mat, x_mat', offspring_hist_mat);
    xlabel('Paternal phenotype'); ylabel('Maternal phenotype');
    zlabel('Offspring phenotype');
    
    
end % if plot
midparent_vec = x_vec(good_bins); offspring_mean_vec = offspring_hist(good_bins); % get output

% Internal function: conditional density of offspring at s given that
% parents phenotype is >= x_mu, for one liability
% s - value of offpring
% h - heritability of each pathway
% x_mu - threshold of parent
%
function ret = conditional_offspring_density_LP(s, h, x_mu)

sigma_t=1;
ret = zeros(1, length(s));
for i=1:length(s)
    ret(i) = quadl(@(t) normpdf(t).*normpdf( (s(i)-h.*t) ./ sqrt(1-h.^2) ), ...
        -10*sigma_t, x_mu) ./ normcdf(x_mu);
end


% Internal function: conditional cumulative of offspring at s given that
% parents phenotype is >= x_mu, for one liability
% s - value of offpring
% h - heritability of each pathway
% x_mu - threshold of parent
%
function ret = conditional_offspring_cumulative_LP(s, h, x_mu)

sigma_t=1;
ret = zeros(1, length(s));
for i=1:length(s)
    ret(i) = quadl(@(t) normpdf(t).*normcdf( (s(i)-h.*t) ./ sqrt(1-h.^2) ), ...
        -10*sigma_t, x_mu) ./ normcdf(x_mu);
end


