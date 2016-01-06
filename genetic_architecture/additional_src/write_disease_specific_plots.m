% Plot all figures for a specific disease
%
% Input:
% disease_specific_plots - flag saying if to plot
% disease_data - more data on .. ???
% data - SNP specific parameters
% data_params - additional parameters
% i - index of trait
% lambda_s_explained_total_percent - how much of sib-risk is explained
% good_familial_inds_vec - do we have epidemiological parameters
% num_shadow_snps_vec - how many shadow snps are there for each snp (cell array with #corrections)
% h_explained_power_adjusted_vec
% lambda_s_var_explained_vec - vector of how much of sib-risk is explained
% root_disease_dir - directory for saving figures
% disease_names - names of traits
% power_correct - matrix of power values (cell array)
% cur_ind_used - index for choosing power
%
% Output:
% var_explained_vec - vector of variance explained by real and shadow snps in decreasing order.
% var_explained_extrapolated_vec - vector of variance explained extrapolated to further loci
% num_loci_needed - estimation of # loci needed for a particular var. explained
% num_snps_inflated - number of snps after power correction (inflation)
% fit_string - best log-log fit to all observed snps
%
function [var_explained_vec var_explained_extrapolated_vec ...
    num_loci_needed num_snps_inflated fit_string ...
    true_beta_mu true_beta_std] = ...
    write_disease_specific_plots(disease_specific_plots, disease_data, data, data_params, i, ...
    lambda_s_explained_total_percent, good_familial_inds_vec, ...
    num_shadow_snps_vec, h_explained_power_adjusted_vec, lambda_s_var_explained_vec, ...
    root_disease_dir, disease_names, power_correct, cur_ind_used)

AssignGeneralConstants;
AssignGeneticArchitectureConstants;
machine = get_machine_type();
fit_string = '';
trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait
num_snps = length(trait_inds);
num_loci_needed = [-1 -1]; % set some default value
%var_explained_vec = []; real_snps_inds_vec = [];
%var_explained_min_vec = []; var_explained_max_vec = [];
%effect_size_vec = [];
var_explained_extrapolated_vec = [];
num_corrections = length(correction_inds);
if(length(cur_ind_used) == 1)
    cur_ind_used = [cur_ind_used ones(1, num_corrections-1)];
end
var_explained_vec = cell(num_corrections,1);
var_explained_min_vec = cell(num_corrections,1);
var_explained_max_vec = cell(num_corrections,1);
real_snps_inds_vec = cell(num_corrections,1);
effect_size_vec = cell(num_corrections,1);
num_snps_inflated_high_power = cell(num_corrections,1);
num_snps_inflated = zeros(num_corrections,1);

for jj = 1:num_corrections % loop on different corrections
    num_snps_inflated_high_power{jj} = find(power_correct{jj}(:,cur_ind_used(jj)) > MIN_POWER); % take only snps with enough power
    num_snps_inflated_high_power{jj} = ...
        sum(num_shadow_snps_vec{jj}(num_snps_inflated_high_power{jj},cur_ind_used(jj))+1)
    
    num_snps_inflated(jj) = num_snps + sum(num_shadow_snps_vec{jj}(:,cur_ind_used(jj)));
    var_explained_vec{jj} = zeros(1, num_snps_inflated(jj));
    var_explained_min_vec{jj} = zeros(1, num_snps_inflated(jj));
    var_explained_max_vec{jj} = zeros(1, num_snps_inflated(jj));
    real_snps_inds_vec{jj} = zeros(1, num_snps_inflated(jj));
    effect_size_vec{jj} = zeros(3, num_snps_inflated(jj)); % allow 3 effect sizes: discovery, replication, combined
    ctr=1;
    for j=1:num_snps % loop on SNPs
        switch correction_mode
            case 'floor'
                cur_num_shadow_snps = floor(num_shadow_snps_vec{jj}(j,cur_ind_used(jj)));
            case 'exact'
                cur_num_shadow_snps = num_shadow_snps_vec{jj}(j,cur_ind_used(jj));
            case 'round'
                cur_num_shadow_snps = round(num_shadow_snps_vec{jj}(j,cur_ind_used(jj)));
        end
        var_explained_vec{jj}(ctr) = ... % [var_explained_vec ... % append real snps
            max(2*epsilon, data_params.snp_h_liab(trait_inds(j))); % ];
        var_explained_vec{jj}(ctr+1:ctr+cur_num_shadow_snps) = ... % [var_explained_vec ...
            repmat(max(epsilon, data_params.snp_h_liab(trait_inds(j))-epsilon), ...
            1, cur_num_shadow_snps); % ]; % append shadow snps - slightly lower (assume an integer number)
        real_snps_inds_vec{jj}(ctr) = 1; % [real_snps_inds_vec 1];
        real_snps_inds_vec{jj}(ctr+1:ctr+cur_num_shadow_snps) = 0;
        
        var_explained_min_vec{jj}(ctr) = ... % [var_explained_min_vec ...
            max(2*epsilon, data_params.snp_h_liab_min(trait_inds(j)));
        var_explained_min_vec{jj}(ctr+1:ctr+cur_num_shadow_snps) = ... % [var_explained_min_vec ...
            repmat(max(epsilon, data_params.snp_h_liab_min(trait_inds(j))-epsilon), ...
            1, cur_num_shadow_snps); % ]; % make shadow snps slightly lower
        var_explained_max_vec{jj}(ctr) = max(2*epsilon, ...
            data_params.snp_h_liab_max(trait_inds(j)));
        var_explained_max_vec{jj}(ctr+1:ctr+cur_num_shadow_snps) = ... %[var_explained_max_vec ...
            repmat(max(epsilon, data_params.snp_h_liab_max(trait_inds(j))-epsilon), ...
            1, cur_num_shadow_snps); % ]; % make shadow snps slightly lower
        
        %     if(j == num_snps_inflated_high_power{jj}) % set number of snps with high confidence
        %         num_snps_inflated_high_power{jj} = length(var_explained_vec);
        %     end
        
        switch disease_data.trait_type{i}
            case {'Binary', 'binary'} % genetic relative risk
                effect_size_vec{jj}(:,ctr) = data.OR_with_MAF(trait_inds(j),:)';
                effect_size_vec{jj}(:,ctr+1:ctr+cur_num_shadow_snps) = ...
                    repmat(data.OR_with_MAF(trait_inds(j),:)-epsilon, ...
                    cur_num_shadow_snps, 1)';
            case {'QTL', 'Quantitative'} % beta. Make sure we've got the correct sign!
                effect_size_vec{jj}(:,ctr) = data.Beta(trait_inds(j),:)';
                effect_size_vec{jj}(:,ctr+1:ctr+cur_num_shadow_snps) = ...
                    repmat(data.Beta(trait_inds(j),:)-epsilon, ...
                    cur_num_shadow_snps, 1)';
        end
        ctr = ctr+cur_num_shadow_snps+1; % num_shadow_snps_vec{jj}(j,cur_ind_used(jj))+1; % add # shadow snps + one real snp
    end % loop on snps
    effect_size_vec{jj} = effect_size_vec{jj}';
    [var_explained_vec{jj} sort_perm] = sort(var_explained_vec{jj}, 'descend');
    var_explained_min_vec{jj} = var_explained_min_vec{jj}(sort_perm);
    var_explained_max_vec{jj} = var_explained_max_vec{jj}(sort_perm);
    [~, effect_sort_perm] = sort(effect_size_vec{jj}(:,2), 'descend'); % sort replication effect sizes: use permutation for variacne explained
    %effect_size_vec_unsorted = effect_size_vec; % keep original ordering of effect sizes
    effect_size_vec{jj} = effect_size_vec{jj}(effect_sort_perm,:); % sort effect sizes
    real_snps_inds_effect_size_vec{jj} = real_snps_inds_vec{jj}(effect_sort_perm);
    shadow_snps_inds_effect_size_vec{jj} = find(~real_snps_inds_effect_size_vec{jj});
    real_snps_inds_effect_size_vec{jj} = find(real_snps_inds_effect_size_vec{jj});
    real_snps_inds_vec{jj} = real_snps_inds_vec{jj}(sort_perm); %        power_vec = max(MIN_POWER, data.Power(trait_inds,:)
    cumsum_var_explained_vec{jj} = cumsum(var_explained_vec{jj});
    shadow_snps_inds_vec{jj} = find(~real_snps_inds_vec{jj});
    real_snps_inds_vec{jj} = find(real_snps_inds_vec{jj});
    real_inds_expand_vec{jj} = zeros(num_snps_inflated(jj),1); % vector of indices that expands the real snps to the shadow ones
    for j=1:num_snps-1
        real_inds_expand_vec{jj}(real_snps_inds_vec{jj}(j):real_snps_inds_vec{jj}(j+1)-1) = j;
    end
    real_inds_expand_vec{jj}(real_snps_inds_vec{jj}(num_snps):end) = num_snps;
end % loop on different power corrections
switch disease_data.trait_type{i} % This may be different from effect_size_vec
    case {'Binary', 'binary'} % genetic relative risk
        observed_effect_size_vec = data.OR_with_MAF(trait_inds,:);
    case {'QTL', 'Quantitative'} % beta. Make sure we've got the correct sign!
        observed_effect_size_vec = data.Beta(trait_inds,:);
end

%num_shadow_snps_vec(j,cur_ind_used)


%num_snps_inflated2 = num_snps + sum(num_shadow_snps_vec(:,cur_ind_used))



switch disease_data.Trait{i}
    case special_traits % {'Height', 'Type 2 diabetes'}
        xxxccc = 12313123
    otherwise
        disease_specific_plots = 0; % plot only for special (to save time ... )
end

true_beta_mu = zeros(num_snps,1); true_beta_std = zeros(num_snps,1); % assign dummy zeros
if(disease_specific_plots)
    x_min = 5; % discard few big effects
    plot_disease_var_vs_sample_size(data, disease_data, data_params, disease_names, ...
        root_disease_dir, i, cur_ind_used, num_snps, ...
        var_explained_vec, real_inds_expand_vec);
    plot_disease_effect_size_hist(disease_data, data_params, disease_names, root_disease_dir, i, ...
        real_snps_inds_vec, shadow_snps_inds_vec, ...
        effect_size_vec, real_snps_inds_effect_size_vec, shadow_snps_inds_effect_size_vec, ...
        num_snps, num_snps_inflated, num_snps_inflated_high_power, cur_ind_used, x_min);
    plot_disease_power_comparison(data, disease_data, data_params, ...
        disease_names, root_disease_dir, i, cur_ind_used);
    
    [true_beta_mu true_beta_std] = ...
        plot_disease_var_by_effect_size(data, disease_data, data_params, disease_names, ...
        root_disease_dir, i, num_snps, num_corrections, cur_ind_used, power_correct);
    
    plot_disease_raw_data(data, disease_data, data_params, disease_names, root_disease_dir, i, ...
        effect_size_vec, observed_effect_size_vec, real_snps_inds_effect_size_vec);
    plot_disease_bar(disease_data, data_params, disease_names, root_disease_dir, i, ...
        good_familial_inds_vec, lambda_s_explained_total_percent);
    plot_disease_cumulative(data, disease_data, data_params, disease_names, root_disease_dir, i, ...
        var_explained_vec, cumsum_var_explained_vec, lambda_s_var_explained_vec, ...
        real_snps_inds_vec, shadow_snps_inds_vec, num_corrections);
    plot_disease_power_sensitivity(disease_data, data_params, disease_names, root_disease_dir, i, ...
        h_explained_power_adjusted_vec);
    area_extrapolated = plot_disease_var_histogram(disease_data, data_params, disease_names, root_disease_dir, i, ...
        real_snps_inds_vec, shadow_snps_inds_vec, ...
        var_explained_vec, cumsum_var_explained_vec, ...
        num_snps, num_snps_inflated, num_snps_inflated_high_power, cur_ind_used, x_min);
    plot_disease_var_cumulative_2x2(disease_data, data_params, disease_names, root_disease_dir, i, ...
        num_snps, num_snps_inflated, num_snps_inflated_high_power, ...
        real_snps_inds_vec, shadow_snps_inds_vec, num_corrections, ...
        var_explained_vec, cumsum_var_explained_vec, cur_ind_used, x_min, area_extrapolated);
    plot_disease_var_individual_2x2(disease_data, data_params, disease_names, root_disease_dir, i, ...
        num_snps, num_snps_inflated, num_snps_inflated_high_power, ...
        real_snps_inds_vec, shadow_snps_inds_vec, num_corrections, num_loci_needed, ...
        var_explained_vec, cumsum_var_explained_vec, ...
        var_explained_min_vec, var_explained_max_vec, cur_ind_used, x_min, area_extrapolated);
end % if disease-specific plots



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal function for specific plots:
function plot_disease_bar(disease_data, data_params, disease_names, root_disease_dir, ...
    i, good_familial_inds_vec, lambda_s_explained_total_percent)
figure; hold on; % Plot a simple bar

switch disease_data.trait_type{i} % show sibling relative risk for binary traits
    case 'Binary'
        total_bar = [1 1];
        familial_bar = [str2double(disease_data.h_familial{i})/100 1]; % str2double(disease_data.lambda_s_familial{i})]
        adjusted_bar = [(data_params.h_liab_power_adjusted(i)/100), ...
            str2double(data_params.lambda_s_power_adjusted_percent{i})/100];
        explained_bar = [data_params.h_liab(i), str2double(lambda_s_explained_total_percent)/100]; %  data_params.lambda_s(i)];
        x_ticks = {'h^2', '\lambda_s'};
    case 'QTL'
        total_bar = 1;
        familial_bar = str2double(disease_data.h_familial{i})/100;
        adjusted_bar = (data_params.h_liab_power_adjusted(i)/100);
        explained_bar = data_params.h_liab(i);
        x_ticks = {'h^2'};
end % switch trait type
bar(total_bar, 'k'); % Plot all disease variance
if(good_familial_inds_vec(i))
    bar(familial_bar, 'r'); % Plot all genetic component
    legend_vec = {'Environmental', 'Genetic-unexplained', 'Explained (power adjusted)', 'Explained (known loci)'};
else
    legend_vec = {'Environmental+Genetic-unexplained', 'Explained (power adjusted)', 'Explained (known loci)'};
end
bar(adjusted_bar, 'b'); % Plot part with power adjustment
bar(explained_bar, 'g'); % Plot part already explained
set(gca, 'XTick', [1 2]); % set labels
set(gca, 'XTickLabel', x_ticks); % set labels
ylabel('h^2'); % Heritability
legend(legend_vec);

switch disease_data.trait_type{i} % show sibling relative risk for binary traits
    case 'Binary'
        lambda_s_ticks = [0 1]; % log(str2double(disease_data.lambda_s_familial{i}))];
        %                lambda_s_ticks = [1:(lambda_s_ticks-1)/10: lambda_s_ticks;
        add_right_yticks(lambda_s_ticks, '\lambda_s'); % right-side scale is for sibling relative risk
        ylabels = get(gca, 'YTickLabel');
        tmp_lambda_familial = str2double(disease_data.lambda_s_familial{i});
        if(isempty(tmp_lambda_familial))
            tmp_lambda_familial = data_params.lambda_s_power_adjusted(i) % temp wrong! (only if no sibling relative risk available!)
        end
        ylabels = num2str(tmp_lambda_familial.^str2double(ylabels), 3);
        set(gca, 'YTickLabel', ylabels);
        %                my_loglog('y'); % set as exponential scale
end
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'genetic_effect_explained_summary'), 'jpg');


% Plot cumulative loci
function plot_disease_cumulative(data, disease_data, data_params, disease_names, root_disease_dir, i, ...
    var_explained_vec, cumsum_var_explained_vec, lambda_s_var_explained_vec, ...
    real_snps_inds_vec, shadow_snps_inds_vec, num_corrections)

trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait

AssignGeneralConstants;
AssignGeneticArchitectureConstants;
figure; hold on; % Plot a cumulative figure with all loci (did we already do that?)
for jj=1:num_corrections
    plot(cumsum_var_explained_vec{jj}(real_snps_inds_vec{jj}), ...
        log(var_explained_vec{jj}(real_snps_inds_vec{jj})), [color_vec(jj) '*']);
    plot(cumsum_var_explained_vec{jj}(shadow_snps_inds_vec{jj}), ...
        log(var_explained_vec{jj}(shadow_snps_inds_vec{jj})), [color_vec(jj) '.']);
end
temp_legend_vec = empty_cell_to_empty_str(cell(2*num_corrections,1));
temp_legend_vec(1:2:end) = data_params.V_corrected_strings{i}(correction_inds);
legend(temp_legend_vec); %     {'Found SNPs', 'Shadow SNPs'});
xlabel('Cumulative variance explained');
ylabel('Individual variance explained (log-scale)');
title('Contribution of discovered loci to total variance of trait');
my_loglog('Y');
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'all_loci_cumulative_effect_size'), 'jpg');

switch disease_data.trait_type{i} % show sibling relative risk for binary traits
    case 'Binary'
        figure; hold on; % Plot var. explained on lambda_s scale vs. liability scale
        plot(100*data_params.snp_h_liab(trait_inds) ./ data.h_familial(trait_inds(1)), ...
            lambda_s_var_explained_vec, '.');
        hold on; plot(lambda_s_var_explained_vec, lambda_s_var_explained_vec, 'r');
        title('Var. explained by each locus on two different scales');
        xlabel('h^2 exp. % (liability-scale)'); ylabel('\lambda_s % explained');
        my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
            'var_explained_two_scales'), 'jpg');
end



% Plot sensitivity of var. explained to computed power
function  plot_disease_power_sensitivity(disease_data, data_params, disease_names, root_disease_dir, i, ...
    h_explained_power_adjusted_vec)


figure; hold on; % Plot sensitivity of var. explained vs. power cutoff
plot(log10(data_params.alpha_vec), 100*sum(h_explained_power_adjusted_vec) ./ ...
    (str2double(disease_data.h_familial{i})/100), '*');
xlabel('p-value cutoff (log_{10})'); ylabel('Var. Explained (power corrected)');
title([str2title(disease_names{i}) '. Sensitivity of Estimated Var. Explained to different Power Corrections']);
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'sensitivity_variance_explained_pval_cutoff'), 'jpg');


% Plot histogram of var. explained of all loci
function area_extrapolated = plot_disease_var_histogram(disease_data, data_params, disease_names, root_disease_dir, i, ...
    real_snps_inds_vec, shadow_snps_inds_vec, ...
    var_explained_vec, cumsum_var_explained_vec, ...
    num_snps, num_snps_inflated, num_snps_inflated_high_power, cur_ind_used, x_min)


figure; hold on; % New figure: histogram of effect sizes (variacne explained)
%         var_explained_vec = sort(100*data_params.snp_h_liab(trait_inds), 'descend');
bar(real_snps_inds_vec{1}, 100*var_explained_vec{1}(real_snps_inds_vec{1}), 'r');  % plot on linear scale (%)
if(~isempty(shadow_snps_inds_vec{1}))
    bar(shadow_snps_inds_vec{1}, 100*var_explained_vec{1}(shadow_snps_inds_vec{1}));  % plot on linear scale (%)
end
area_extrapolated = 0; % default: no extrapolation
if((num_snps > 10) && (num_snps_inflated_high_power{1} > x_min)) % no point in fitting too few points
    exp_fit = polyfit(x_min:num_snps_inflated_high_power{1}, ...
        log(var_explained_vec{1}(x_min:num_snps_inflated_high_power{1})), 1)  %fit an exponential
    %            exp_fit2 = fit(vec2column(x_min:num_snps_inflated), vec2column(var_explained_vec(x_min:end)), 'exp1')
    powerlaw_density_fit = polyfit(log(x_min:num_snps_inflated_high_power{1})', ...
        log(var_explained_vec{1}(x_min:num_snps_inflated_high_power{1})'), 1)
    %            powerlaw_density_fit = fit((x_min:num_snps_inflated)', var_explained_vec(x_min:num_snps_inflated)', 'power2')
    
    plot(x_min:num_snps_inflated(1)*1.1, 100 * exp(exp_fit(2)) .* exp([x_min:num_snps_inflated(1)*1.1].*(exp_fit(1))), ...
        'g', 'linewidth', 2); % plot exponential
    %            plot(x_min:num_snps_inflated*1.1, 100 * exp_fit2.a *  exp(exp_fit2.b .* [x_min:num_snps_inflated*1.1]), ...
    %                'c', 'linewidth', 2); % plot exponential
    %            plot(x_min:num_snps_inflated*1.1, 100 * (powerlaw_density_fit.a * ...
    %                [x_min:num_snps_inflated*1.1].^powerlaw_density_fit.b + powerlaw_density_fit.c), ...
    %                'm', 'linewidth', 2); % plot power-law
    plot(x_min:num_snps_inflated(1)*1.1, 100 * exp(powerlaw_density_fit(2)) * ...
        [x_min:num_snps_inflated(1)*1.1].^powerlaw_density_fit(1), ...
        'm', 'linewidth', 2); % plot power-law
    area_extrapolated = 100*(-exp(exp_fit(2)) / exp_fit(1)) * ...
        exp(exp_fit(1) * num_snps_inflated(1)); % add extrapolation
end
title(['Var. Explained (Sorted). Total: ' num2str(data_params.h_liab(i)) ...
    '% (observed), ' num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i)) ...
    '% (inflated), ' ...
    num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i)+area_extrapolated) ...
    '% (extrapolated)'], 'fontsize',12);
xlabel('rank'); ylabel('Var. Explained (%)');
legend({'found SNPs', 'shadow SNPs'}); % put on bottom right corner
yl = ylim;
if(yl(2) > 0) % avoid problems ...
    ylim([0,  min(2, yl(2))]); % don't show first largest effects to not skew effects
end
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'all_loci_variance_explained'), 'jpg');

% plot cumulative var. explained vs. ... on linear and log scale
function    plot_disease_var_cumulative_2x2(disease_data, data_params, disease_names, root_disease_dir, i, ...
    num_snps, num_snps_inflated, num_snps_inflated_high_power, ...
    real_snps_inds_vec, shadow_snps_inds_vec, num_corrections, ...
    var_explained_vec, cumsum_var_explained_vec, cur_ind_used, x_min, area_extrapolated)

AssignGeneralConstants;
machine = get_machine_type();
switch machine
    case UNIX
        bmp_fig_format = {'jpg'};
    case PC
        bmp_fig_format = {'bmp', 'jpg'};
end

do_exp_fit = 0; % dont fit exponentials (just power-law)


full_figure; % Make a big 2X2 plot of cumulative variance
lots_snps_flag = double(num_snps > 10); % run only if many SNPs
for j=0:lots_snps_flag
    for k=1:2 % do on linear and log scale as well
        subplot(2,2,2*j+k); hold on;  %       figure; hold on; % New figure: cumulative plots of effect sizes (variacne explained)
        if(k == 1)
            x_vec = 1:num_snps_inflated(1);
        else
            x_vec = log10(1:num_snps_inflated(1));
        end
        plot(x_vec(real_snps_inds_vec{1}), ...
            100*cumsum_var_explained_vec{1}(real_snps_inds_vec{1}), 'ro');  % plot on linear scale (%)
        if(~isempty(shadow_snps_inds_vec))
            cumsum_real_var_explained_vec{1} = cumsum(var_explained_vec{1}(real_snps_inds_vec{1})); % sum up ONLY real SNPs
            plot(x_vec(shadow_snps_inds_vec{1}), 100*cumsum_var_explained_vec{1}(shadow_snps_inds_vec{1}), '.');  % plot on linear scale (%)
            plot(x_vec(1:num_snps), 100*cumsum_real_var_explained_vec{1}, 'r.');  % plot again only the real SNPs.
        end
        if((num_snps > 10) && (num_snps_inflated_high_power{1} > x_min))% no point in fitting too few points
            % Fit an power-law curve to cumulative points
            powerlaw_fit = polyfit(log(x_min:num_snps_inflated_high_power{1})', ...
                log(var_explained_vec{1}(x_min:num_snps_inflated_high_power{1})'), 1)
            ff = fittype([num2str(exp(powerlaw_fit(2))/(powerlaw_fit(1)+1)) ...
                '*x^' num2str(powerlaw_fit(1)+1) '+c']); % try different fitting options
            powerlaw_fit2 = fit((x_min:num_snps_inflated(1))', ...
                cumsum_var_explained_vec{1}(x_min:num_snps_inflated(1))', ff)
            if(do_exp_fit)
                exponential_fit3 =  polyfit((x_min:num_snps_inflated_high_power{1})', ...
                    log(var_explained_vec{1}(x_min:num_snps_inflated_high_power{1}))', 1)
                exponential_fit.a = exp(exponential_fit3(2));
                exponential_fit.b = exponential_fit3(1);
                ff = fittype([num2str(exponential_fit.a/exponential_fit.b) ...
                    '*exp(x*' num2str(exponential_fit.b) ')+c']); % try different fitting options
                exponential_fit2 = fit((x_min:num_snps_inflated(1))', ...
                    cumsum_var_explained_vec{1}(x_min:num_snps_inflated(1))', ff)
            end
            if(j == 0)
                x_vec = x_min:num_snps_inflated(1); % Make a large plot    %   x_min:num_snps_inflated(1)*1.1
            else
                x_vec = 1:5000;
                half_heritability = str2double(disease_data.h_familial{i}) / 200;
                num_loci_needed(1) = ((half_heritability - powerlaw_fit2.c) * (1+powerlaw_fit(1)) / exp(powerlaw_fit(2))) ^ (1 / (1+powerlaw_fit(1))); % num needed by exponential distribution
                if(do_exp_fit)
                    num_loci_needed(2) = log((half_heritability - exponential_fit2.c) / (exponential_fit.a/exponential_fit.b)) / exponential_fit.b; % num needed by exponential distribution
                end
            end
            y_power_vec = max(0, 100*(powerlaw_fit2.c + exp(powerlaw_fit(2))/(powerlaw_fit(1)+1).*x_vec.^(powerlaw_fit(1)+1)));
            if(do_exp_fit)
                y_exp_vec = max(0, 100*(exponential_fit2.c + (exponential_fit.a/exponential_fit.b).*exp(exponential_fit.b .* (x_vec))));
            end
            if(j <= 1)
                y_power_vec = min(y_power_vec, str2double(disease_data.h_familial{i}));
                var_explained_extrapolated_vec = 100 .* y_power_vec ./ ...
                    str2double(disease_data.h_familial{i});
                if(do_exp_fit)
                    y_exp_vec = min(y_exp_vec, str2double(disease_data.h_familial{i}));
                end
            end
            if(k == 2)
                x_vec = log10(x_vec);
            end
            plot(x_vec,  y_power_vec, 'g', 'linewidth', 2); % plot power-law
            if(do_exp_fit)
                plot(x_vec,  y_exp_vec, 'm', 'linewidth', 2); % plot exponential
            end
            if((k == 2) && (j == 1))
                if(do_exp_fit)
                    legend({'found SNPs', 'shadow SNPs', '', 'pow-law-fit', 'exp-fit'}, 4, 'fontsize', 8);
                else
                    legend({'found SNPs', 'shadow SNPs', '', 'pow-law-fit'}, 4, 'fontsize', 8);
                end
            end
        else
            legend({'found SNPs', 'shadow SNPs'}, 4);
        end % if num_snps > 10
        log_str = repmat(' (log_{10}) ', 1, k-1);
        if(j == 0)
            cum_fig_name = 'all_loci_cumulative_variance_explained';
            if(k == 1)
                title(['Cum. Var. explained. Tot.: ' num2str(100*data_params.h_liab(i),3) ...
                    '% (obs.), ' num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i),3) ...
                    '% (inflated), ' ...
                    num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i)+ ...
                    area_extrapolated,3) '% (extrap.)'], 'fontsize',12);
            else % Plot function fitted on the right side
                if( (num_snps > 10) && (num_snps_inflated_high_power{1}>x_min) )
                    fit_string = ['V<sub>n</sub> = ' num2str(exp(powerlaw_fit(2)),3) ...
                        '*n<sup>' num2str(powerlaw_fit(1),3) '</sup>']; % html format
                    title(['V(n) = ' num2str(exp(powerlaw_fit(2)),3) '*n^{' ...
                        num2str(powerlaw_fit(1),3) '}. ' ...
                        'Cum. V(n) = ' num2str(powerlaw_fit2.c,3) ' + ' ...
                        num2str(exp(powerlaw_fit(2))/(powerlaw_fit(1)+1),3) ...
                        '*n^{' num2str((powerlaw_fit(1)+1),3) '}'], 'fontsize',12);
                else
                    fit_string = ''; % no fitting ...
                    title('no fitting (too few loci ..)');
                end
            end
        else
            %            cum_fig_name = 'all_loci_cumulative_variance_explained_extrapolation';
            if(k == 1)
                title(['Cum. h^2 Explained. ' ...
                    num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i)+area_extrapolated, 3) '% (extrapolated). ' ...
                    'N=' num2str(round(num_loci_needed(1))) ' loci to reach 50% h^2'], 'fontsize', 12); % report #loci needed based on power-law fit
            end
        end
        xlabel(['rank' log_str]);
        if(k == 1)
            ylabel('Cumulative Var. Explained (%)');
        end
        y_ticks = get(gca, 'YTick'); % y_ticks = [y_ticks 2*y_ticks(end) - y_ticks(end-1)]; % add last element
        y_ticks = y_ticks*100 ./ str2double(disease_data.h_familial{i});
        if(k == 2)
            y_right_label = 'Cumulative h^2 Explained (%)';
        else
            y_right_label = '';
        end
        add_right_yticks([y_ticks(1) y_ticks(end)], y_right_label); % set scale to heritability (need to round it)
    end % loop on linear vs. log scale
end % loop on two extrapolation plots
mtit(str2title(disease_names{i}), 'fontsize',16,'color','r'); % a major title for all sub-plots
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    cum_fig_name), bmp_fig_format); % jpg screws right hand size axis



% Plot individual var. explained (vs. cumulative?) on log and linear scales
function       plot_disease_var_individual_2x2(disease_data, data_params, disease_names, root_disease_dir, i, ...
    num_snps, num_snps_inflated, num_snps_inflated_high_power, ...
    real_snps_inds_vec, shadow_snps_inds_vec, num_corrections, num_loci_needed, ...
    var_explained_vec, cumsum_var_explained_vec, ...
    var_explained_min_vec, var_explained_max_vec, cur_ind_used, x_min, area_extrapolated)


AssignGeneralConstants;
machine = get_machine_type();
switch machine
    case UNIX
        bmp_fig_format = {'jpg'};
    case PC
        bmp_fig_format = {'bmp', 'jpg'};
end


full_figure; % Make a big 2X2 plot of cumulative and individual variance on log scale
for j=1:2 % here just do cumulative vs. marginal variance
    for k=1:2 % do x axis on log scale as well
        if((k == 2) && (j == 1)) % make two plots
            num_plots = 2;
        else
            num_plots = 1;
        end
        for cur_plot=1:num_plots
            if(cur_plot == 1)
                subplot(2,2,2*(j-1)+k); hold on;  %       figure; hold on; % New figure: cumulative plots of effect sizes (variacne explained)
            else
                full_figure;
            end
            if(k == 1)
                x_vec = 1:num_snps_inflated(1);
            else
                x_vec = log10(1:num_snps_inflated(1));
            end
            if(j == 1) % individual effects
                y_vec = log10(var_explained_vec{1});
            else % cumulative effects
                y_vec = log10(cumsum_var_explained_vec{1});
            end
            plot(x_vec(real_snps_inds_vec{1}), y_vec(real_snps_inds_vec{1}), 'r.');  % plot on linear scale (%)
            if(j == 1) % New: Plot also confidence intervals
                errorbar(x_vec(real_snps_inds_vec{1}), y_vec(real_snps_inds_vec{1}), ...
                    y_vec(real_snps_inds_vec{1}) - log10(var_explained_min_vec{1}(real_snps_inds_vec{1})), ...
                    log10(var_explained_max_vec{1}(real_snps_inds_vec{1})) - y_vec(real_snps_inds_vec{1}));
            end
            if(~isempty(shadow_snps_inds_vec{1}))
                %                cumsum_real_var_explained_vec = cumsum(var_explained_vec(real_snps_inds_vec)); % sum up ONLY real SNPs
                plot(x_vec(shadow_snps_inds_vec{1}), y_vec(shadow_snps_inds_vec{1}), '.');  % plot on logarithmic scale (%)
                %                plot(x_vec(1:num_snps), 100*cumsum_real_var_explained_vec, 'ro');  % plot again only the real SNPs.
            end
            if( (num_snps > 10) && (num_snps_inflated_high_power{1}>x_min) )% no point in fitting too few points
                % Fit an power-law curve to cumulative points
                powerlaw_fit = polyfit(log(x_min:num_snps_inflated_high_power{1})', ...
                    log(var_explained_vec{1}(x_min:num_snps_inflated_high_power{1})'), 1)
                ff = fittype([num2str(exp(powerlaw_fit(2))/(powerlaw_fit(1)+1)) '*x^' num2str(powerlaw_fit(1)+1) '+c']); % try different fitting options
                powerlaw_fit2 = fit((x_min:num_snps_inflated(1))', ...
                    cumsum_var_explained_vec{1}(x_min:num_snps_inflated(1))', ff)
                x_vec = x_min:num_snps_inflated(1); % Make a large plot    %   x_min:num_snps_inflated(1)*1.1
                
                if(j == 1) % plot marginal var. explained
                    y_power_vec = max(0, (exp(powerlaw_fit(2)).*x_vec.^(powerlaw_fit(1))));
                else % plot cumulative var. explained
                    y_power_vec = max(0, (powerlaw_fit2.c + exp(powerlaw_fit(2))/(powerlaw_fit(1)+1).*x_vec.^(powerlaw_fit(1)+1)));
                end
                y_power_vec = min(y_power_vec, str2double(disease_data.h_familial{i}));
                if(k == 2)
                    x_vec = log10(x_vec);
                end
                plot(x_vec, log10(y_power_vec), 'g', 'linewidth', 2); % plot power-law
                if((k == 2) && (j == 2))
                    legend({'found SNPs', 'shadow SNPs', 'pow-law-fit'}, 4, 'fontsize', 8);
                end
            else
                legend({'found SNPs', 'shadow SNPs'}, 4);
            end % if num_snps > 10
            log_str = repmat(' (log_{10}) ', 1, k-1);
            if(j == 1)
                cum_fig_name = 'all_loci_variance_explained_log_scale';
                if(k == 1)
                    title(['Cum. Var. explained. Tot.: ' num2str(100*data_params.h_liab(i),3) ...
                        '% (obs.), ' num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i),3) '% (inflated), ' ...
                        num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i)+area_extrapolated,3) '% (extrap.)'], 'fontsize',12);
                else % Plot function fitted on the right side
                    if((num_snps > 10) && (num_snps_inflated_high_power{1}>x_min))
                        title(['V(n) = ' num2str(exp(powerlaw_fit(2)),3) '*n^{' num2str(powerlaw_fit(1),3) '}. ' ...
                            'Cum. V(n) = ' num2str(powerlaw_fit2.c,3) ' + ' ...
                            num2str(exp(powerlaw_fit(2))/(powerlaw_fit(1)+1),3) '*n^{' num2str((powerlaw_fit(1)+1),3) '}'], 'fontsize',12);
                    else
                        title('No fitting (too few SNPs ...)');
                    end
                end % if k == 1
            else
                %            cum_fig_name = 'all_loci_cumulative_variance_explained_extrapolation';
                if(k == 1)
                    title(['Cum. h^2 Explained. ' ...
                        num2str(data_params.h_liab_power_adjusted(cur_ind_used(1),i)+area_extrapolated, 3) '% (extrapolated). ' ...
                        'N=' num2str(round(num_loci_needed(1))) ' loci to reach 50% h^2'], 'fontsize',12); % report #loci needed based on power-law fit
                end
            end
            if(j == 1)
                ylabel('Var. Explained log_{10}.');
            else
                ylabel('Cumulative Var. Explained log_{10}.');
            end
            y_ticks = get(gca, 'YTick'); % y_ticks = [y_ticks 2*y_ticks(end) - y_ticks(end-1)]; % add last element
            y_ticks = y_ticks*100 ./ str2double(disease_data.h_familial{i});
            if(k == 2)
                if(j == 1)
                    y_right_label = 'h^2 Explained log_{10}';
                else
                    y_right_label = 'Cumulative h^2 Explained log_{10}';
                end
            else
                y_right_label = '';
            end
            add_right_yticks([y_ticks(1) y_ticks(end)], y_right_label); % set scale to heritability (need to round it)
            
            if(cur_plot == 2)
                xlabel([str2title(disease_names{i}) ', rank' log_str]);
                my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
                    [cum_fig_name '_loglog']), bmp_fig_format); % jpg screws right hand size axis
                close; % close current figure
            else
                xlabel(['rank' log_str]);
            end
        end % loop on cur_plot
    end % loop on linear vs. log scale
end % loop on two extrapolation plots

mtit(str2title(disease_names{i}), 'fontsize',16,'color','r'); % a major title for all subplots
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    cum_fig_name), bmp_fig_format); % jpg screws right hand size axis


% Plot var. explained by discovered loci vs. study sample size
function   plot_disease_var_vs_sample_size(data, disease_data, data_params, disease_names, ...
    root_disease_dir, i, cur_ind_used, num_snps, ...
    var_explained_vec, real_inds_expand_vec)

AssignGeneralConstants;
machine = get_machine_type();
switch machine
    case UNIX
        bmp_fig_format = {'jpg'};
    case PC
        bmp_fig_format = {'bmp', 'jpg'};
end

trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait
full_figure(); % hold on; % New figure: plot var. explained vs. sample size
n_sample_size = round(10.^[0.1:0.1:0.6]./2).*2; %% round(10.^[0.1:0.1:6]./2).*2; % 1000:1000:100000;
num_n_samples = length(n_sample_size); % vector of different sample sizes
alpha = data_params.alpha_vec(cur_ind_used(1)); %   5*10^(-8);
p_z_x_marginal = genetic_relative_risk_to_p_z_x_marginal(data.RAF(trait_inds), ...
    data.OR(trait_inds,2), data.Prevalence(trait_inds)); % get 4XN matrix using REPLICATION effect sizes !!!
%    p_z_x_marginal = p_z_x_marginal(real_inds_expand_vec,:);
%    sum_discovered_var_explained = zeros(num_n_samples, 1);
power_sample_vec = zeros(num_n_samples, num_snps);
for j=1:num_snps % compute power again for various sample sizes.
    switch data.trait_type{trait_inds(1)}
        case 'Binary' % test for a binary SNP
            power_sample_vec(:,j) = ...
                compute_association_power(p_z_x_marginal(j,:), ...
                n_sample_size, [], alpha, ...
                [], 'armitage', 'chi-square-analytic');
        case {'Quantitative', 'QTL'} % test for QTL
            power_sample_vec(:,j) = ...
                compute_association_power(p_z_x_marginal(j,:), ...
                n_sample_size./2, n_sample_size./2, alpha, ...
                [], 'single-locus', 'chi-square-QTL-analytic', 'population');
    end
end
power_sample_vec = power_sample_vec(:,real_inds_expand_vec{1}); % expand to include shadow snps
sum_discovered_var_explained = sum(power_sample_vec .* ...
    repmat(var_explained_vec{1}, num_n_samples,1),2);
num_plots = 1;
for j=1:num_plots
    switch j
        case num_plots-2
            x_vec = n_sample_size;
            x_str = '';
        case num_plots-1
            x_vec = sqrt(n_sample_size);
            x_str = ' (sqrt) ';
        case num_plots
            x_vec = log10(n_sample_size);
            x_str = ' (log_{10}-scale) ';
    end
    subplot(num_plots,1,j); hold on;
    plot(x_vec, 100*sum_discovered_var_explained, '.'); % plot points
    plot(x_vec, 100*sum_discovered_var_explained, 'r'); % interpolation
    line(log10([data_params.effective_sample_size(trait_inds(1)) data_params.effective_sample_size(trait_inds(1))]), ...
        [0 sum(data_params.snp_h_liab(trait_inds))*100], ...
        'color', 'k');
    line(log10([1 data_params.effective_sample_size(trait_inds(1))]), ...
        [sum(data_params.snp_h_liab(trait_inds))*100 sum(data_params.snp_h_liab(trait_inds))*100], ...
        'color', 'k');
    
    xlabel(['N samples ' x_str]); ylabel('Cum. Var. Explained (%)');
    my_title(['Var explained vs. sample size for ' disease_names{i}]);
    x_ticks = get(gca, 'XTick');
    x_tick_labels = num2str_cell(num2cell(10.^x_ticks))
    set(gca, 'XTickLabel', x_tick_labels)
    y_ticks = get(gca, 'YTick'); % y_ticks = [y_ticks 2*y_ticks(end) - y_ticks(end-1)]; % add last element
    y_ticks = y_ticks*100 ./ str2double(disease_data.h_familial{i});
    add_right_yticks([y_ticks(1) y_ticks(end)], 'Cumulative h^2 Explained (%)'); % set scale to heritability (need to round it)
end
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'var_explained_vs_sample_size'), bmp_fig_format); % jpg screws right hand size axis



% Plot histogram of effect sizes (beta or GRR)
function    plot_disease_effect_size_hist(disease_data, data_params, disease_names, root_disease_dir, i, ...
    real_snps_inds_vec, shadow_snps_inds_vec, ...
    effect_size_vec, real_snps_inds_effect_size_vec, shadow_snps_inds_effect_size_vec, ...
    num_snps, num_snps_inflated, num_snps_inflated_high_power, cur_ind_used, x_min)


figure; hold on; % New figure: histogram of effect sizes (genetic-relative-risk or beta)
switch disease_data.trait_type{i} % bar plots. like before but different X axis
    case 'Binary'
        % Get also shadow SNPs
        bar(real_snps_inds_effect_size_vec{1}, ...
            log(effect_size_vec{1}(real_snps_inds_effect_size_vec{1},2)), 'r');
        if(~isempty(shadow_snps_inds_effect_size_vec{1}))
            bar(shadow_snps_inds_effect_size_vec{1}, ...
                log(effect_size_vec{1}(shadow_snps_inds_effect_size_vec{1},2)));
        end
        
        title('Genetic Relative Risk (Sorted)');
        xlabel('rank'); ylabel('GRR (log-scale)');
        legend('found SNPs', 'shadow SNPs', 'extrapolation');
        if((num_snps > 10) && (num_snps_inflated_high_power{1}>x_min))% no point in fitting too few points
            x_min = 5; % discard few big effects
            exp_fit = polyfit(x_min:num_snps_inflated_high_power{1}, ...
                vec2row(log(log(effect_size_vec{1}(x_min:num_snps_inflated_high_power{1},2)))), 1)  %fit an exponential
            plot(x_min:(num_snps_inflated(1)*1.1), exp(exp_fit(2)) .* exp([x_min:(num_snps_inflated(1)*1.1)].*(exp_fit(1))), ...
                'g', 'linewidth', 2); % plot exponential
            %                    figure; hold on; plot(x_min:num_snps_inflated(1), log(effect_size_vec(x_min:end)), '.')
            %                    plot(x_min:num_snps_inflated(1), exp_fit(2) + exp_fit(1) .* (x_min:num_snps_inflated(1)), ...
            %                        'g', 'linewidth', 2); % plot exponential
            %                    figure; hold on; plot(x_min:num_snps_inflated(1), effect_size_vec(x_min:end), '.')
            %                     [alpha_effect, xmin_effect, L_effect] = plfit(effect_size_vec); % fit a power-low
            %                     plot(1:num_snps_inflated(1), [1:num_snps_inflated(1)].^(-alpha_effect), 'g'); % plot power-law
        end
        my_loglog('y'); % plot on log-scale
        
        my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
            'all_loci_genetic_relative_risk'), 'jpg');
    case {'QTL', 'Quantitative'}
        bar(real_snps_inds_effect_size_vec{1}, ...
            effect_size_vec{1}(real_snps_inds_effect_size_vec{1},2), 'r'); % my_loglog('y'); % plot on absolute scale
        if(~isempty(shadow_snps_inds_effect_size_vec{1}))
            bar(shadow_snps_inds_effect_size_vec{1}, ...
                effect_size_vec{1}(shadow_snps_inds_effect_size_vec{1},2)); % my_loglog('y'); % plot on absolute scale
        end
        title('QTL Effect size \beta (Sorted)');
        xlabel('rank'); ylabel('\beta');
        legend('found SNPs', 'shadow SNPs', 'extrapolation');
        x_min = 5; % discard few big effects
        if((num_snps > 10) &&  (num_snps_inflated_high_power{1} > x_min))% no point in fitting too few points
            exp_fit = polyfit(x_min:num_snps_inflated_high_power{1}, ...
                vec2row(log(effect_size_vec{1}(x_min:num_snps_inflated_high_power{1},2))), 1)  %fit an exponential
            plot(x_min:num_snps_inflated(1)*1.1, exp(exp_fit(2)) .* ...
                exp([x_min:num_snps_inflated(1)*1.1].*(exp_fit(1))), ...
                'g', 'linewidth', 2); % plot exponential
            %                     [alpha_effect, xmin_effect, L_effect] = plfit(effect_size_vec); % fit a power-low
            %                     plot(1:num_snps_inflated(1), [1:num_snps_inflated(1)].^(-alpha_effect), 'g'); % plot power-law
        end
        my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
            'all_loci_beta_effect_size'), 'jpg');
end % switch disease trait type

% Compare empirical to theoretical power calculations
function    plot_disease_power_comparison(data, disease_data, data_params, ...
    disease_names, root_disease_dir, i, cur_ind_used)

trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    switch disease_data.trait_type{i} % bar plots. like before but different X axis
%        case 'Binary'
figure; hold on; % Plot two power calculations figure
plot(data.Power(trait_inds,cur_ind_used(1)), ...
    data_params.Power_empirical(trait_inds,cur_ind_used(1)), '.');
plot(0:0.1:1, 0:0.1:1, 'r');
my_title(['Comparison of theoretical vs. empirical powers for ' disease_names{i}]);
xlabel('Theoretical power'); ylabel('Empirical power');
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'power_empirical_vs_theoretical'), 'jpg');
%    end


% Plot raw effect sizes vs. MAF (discovery and replicated)
function plot_disease_raw_data(data, disease_data, data_params, disease_names, root_disease_dir, i, ...
    effect_size_vec, observed_effect_size_vec, real_snps_inds_effect_size_vec)

full_figure; % New: plot 'raw' effect sizes, with power curve
%    N=10^3; x_grid = 0.5.*(1:N)./N; % set effective population size and MAF grid (small N for faster running time)
trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait

switch data.trait_type{trait_inds(1)} % convension: always take MAF !
    case {'binary', 'Binary'}
        %            [x_vec true_grr_vec] = flip_allele(x_vec, true_grr_vec, 'binary', 1);
        %            true_beta_vec = [observed_effect_size_vec(:,1) repmat(data.Prevalence(trait_inds(1)), num_snps, 1)];
        test_type = 'armitage'; test_stat = 'chi-square';
    case {'QTL', 'Quantitative'}
        %           [x_vec true_beta_vec] = flip_allele(x_vec, true_beta_vec, 'QTL', 1);
        test_type = 'single-locus'; test_stat = 'chi-square-QTL';
        %            num_controls = [];
end
if(~isfield(data_params, 'beta_discovery_grid'))
    [beta_discovery_grid beta_discovery_inv_grid] = ...
        compute_discovery_boundary(data_params.x_grid, data.Prevalence(trait_inds(1)), ...
        data.discovery_num_cases(trait_inds(1)), ...
        data.discovery_num_controls(trait_inds(1)), alpha, ...
        test_type, test_stat); % Determine discovery boundry (takes some time)
else
    beta_discovery_grid = data_params.beta_discovery_grid{i};
    beta_discovery_inv_grid = data_params.beta_discovery_inv_grid{i};
end

legend_vec = [];
if(isfield(data, 'OR')) %%%% this is the DISCOVERY effect
    plot(data.MAF(trait_inds), ...
        observed_effect_size_vec(:,1), 'ro'); % 1 is discovery effect
    legend_vec = [legend_vec {'Found SNPs (discovery effect size)'}];
end
if(isfield(data, 'OR')) % this is the REPLICATION effect size
    plot(data.MAF(trait_inds), ...
        observed_effect_size_vec(:,2), 'r+'); % 2 is replication effect
    legend_vec = [legend_vec {'Found SNPs (replication effect size)'}];
end

plot(data_params.x_grid, [beta_discovery_grid beta_discovery_inv_grid], 'g--', 'linewidth', 3); % plot discovery boundary
title(['Effect Size vs. Allele. Freq. All Loci, ' data.Trait{trait_inds(1)}]);
xlabel('MAF');
legend_vec = [legend_vec {'discovery-boundary'}];
legend(legend_vec);
% Don't show bins
%     for j=1:length(bin_boundaries)
%         line([0 0.5], [bin_boundaries(j) bin_boundaries(j)], ...
%             'color', 'k', 'linestyle', '--');
%     end
switch data.trait_type{trait_inds(1)}
    case {'Binary', 'binary'}
        plot(data_params.x_grid, repmat(1, length(data_params.x_grid), 1), 'k:'); % plot line GRR=1
        y_min = 0; effect_str = 'GRR';
    case {'QTL', 'Quantitative'}
        y_min = 0; effect_str = '\beta';
end
ylabel(effect_str);
ylim([y_min max(y_min+0.1, min(max(max(abs(effect_size_vec{1}(real_snps_inds_effect_size_vec{1},1)), ...
    abs(effect_size_vec{1}(real_snps_inds_effect_size_vec{1},2))))*1.1))]); % (50) make scale smaller
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'effect_size_vs_MAF'), 'jpg');




% Plot histogram of variance explained for different effect sizes
function [true_beta_mu true_beta_std] = ...
    plot_disease_var_by_effect_size(data, disease_data, data_params, disease_names, ...
    root_disease_dir, i, num_snps, num_corrections, cur_ind_used, power_correct)

AssignGeneralConstants;
AssignGeneticArchitectureConstants;
alpha = data_params.alpha_vec(cur_ind_used(1)); %   5*10^(-8);
true_beta_mu = zeros(num_snps,1);
true_beta_std = zeros(num_snps,1);


trait_inds = strmatch(data_params.trait_name{i}, data.Trait, 'exact'); % get all indices of current trait
switch data.trait_type{trait_inds(1)}
    case {'Binary', 'binary'}
        effect_str = 'GRR';
    case {'QTL', 'Quantitative'}
        effect_str = '\beta';
end

full_figure;  % New: figure of var. explained dependency on beta
switch data.trait_type{trait_inds(1)}
    case {'QTL', 'Quantitative'}
        beta_grid = -1:0.001:1; % all possible betas
        test_type = 'single-locus';
        test_stat = 'chi-square-QTL';
        
    case {'Binary', 'binary'}
        beta_grid = 0:0.001:10; % all possible GRR's
        test_type = 'armitage';
        test_stat = 'chi-square';
        
end
%    num_corrections = 2; % so far try only one correction
all_beta_estimate_corrected = zeros(num_corrections+1, length(beta_grid));
data_params.Correction = zeros(data_params.num_snps, num_corrections+1);
data_params.Correction(trait_inds,1) = 1; % no correction at all
for jj=1:num_corrections
    data_params.Correction(trait_inds,jj+1) = 1 ./ power_correct{jj}(:,cur_ind_used(jj));
    %        data_params.Power(trait_inds,cur_ind_used(1)); % Temp!!! take the same correction!!!
end
for j=1:num_snps % loop and add contribution to variance explained
    cur_x = data.RAF(trait_inds(j));
    switch data.trait_type{trait_inds(1)} % flip allele
        case {'QTL', 'Quantitative'}
            cur_replication_beta = data.Beta(trait_inds(j),2); % just copy effect size
        case {'Binary', 'binary'}
            cur_replication_beta = [data.OR(trait_inds(j),2) data.Prevalence(trait_inds(j))]; % just copy effect size
            if(cur_replication_beta(1) < 1)
                cur_replication_beta(1) = 1 / cur_replication_beta(1);
                cur_x = 1-cur_x;
            end
    end
    [cur_true_beta_estimate true_beta_mu(j) true_beta_std(j)] = ... % estimation of true beta (indep. of correction used)
        observed_effect_size_dist(cur_replication_beta, ...
        cur_x, test_stat, alpha, beta_grid, ...
        data.replication_num_cases(trait_inds(j)), ...
        data.replication_num_controls(trait_inds(j)), 'replication');
    
    for jj=1:num_corrections+1 % look how much is corrected for each correction type
        all_beta_estimate_corrected(jj,:) = ... % still need to fill in data_params.Correction
            all_beta_estimate_corrected(jj,:) + ...
            cur_true_beta_estimate .* data_params.snp_h_liab(trait_inds(j)) .* ...
            data_params.Correction(trait_inds(j),jj); % corrected variance (pointwise)
    end
end
switch data.trait_type{trait_inds(1)} % collapse MAF [0,0.5] and [0.5,1]
    case 'QTL' % flip beta to collapse positive and negative together
        mid_ind = (length(beta_grid)-1)/2;
        beta_grid = beta_grid(mid_ind+2:end);
        all_beta_estimate_corrected = all_beta_estimate_corrected(:,mid_ind+2:end) + ...
            all_beta_estimate_corrected(:,mid_ind:-1:1);
    case {'Binary', 'binary'}
        mid_ind = find(beta_grid < 1, 1, 'last');
        %                                beta_inv_grid = min(100, 1 ./ beta_grid(1:mid_ind)); % don't get too high GRR
        beta_grid = beta_grid(mid_ind+2:end);
        all_beta_estimate_corrected = ...
            all_beta_estimate_corrected(:,mid_ind+2:end);
        %                                for jj=1:num_corrections
        %                                    [temp_beta_grid temp_all_beta_estimate_corrected(jj,:)] = ...
        %                                    sum_hist(beta_grid, all_beta_estimate_corrected(jj,mid_ind+2:end), ...
        %                                    beta_inv_grid(1:mid_ind), ...
        %                                    all_beta_estimate_corrected(jj,1:mid_ind), [], 0);
        %                                end
        %                                beta_grid = temp_beta_grid;
        %                                all_beta_estimate_corrected = temp_all_beta_estimate_corrected;
end
show_inds = 1:num_corrections+1; % temp: show all corrections
line_width_vec = ones(length(show_inds), 1); line_width_vec(1) = 3; % make true effect size thicker
for jj=1:length(show_inds) % num_corrections. We should include the TWO corrections AND original data!!!
    %        if(show_inds(jj) <= num_corrections)
    plot(beta_grid, all_beta_estimate_corrected(show_inds(jj),:), ...
        color_vec(jj), 'linewidth', line_width_vec(jj));
    %                              total_V_explained_integral(jj) = ...
    %                                  integral_hist(beta_grid, ...
    %                                  all_beta_estimate_corrected(show_inds(jj),:));
end % loop on difference corrections
title(['Var. explained as function of effect size, ' str2title(disease_names{i}) ]);
xlabel(effect_str); ylabel('Var. Explained'); % should be '|\beta|'
legend_vec = ['Observed' data_params.V_corrected_strings{i}(correction_inds)]; %  {'Observed', 'Corrected '};
cur_legend_vec = cell(length(show_inds),1);
total_V_explained_integral = zeros(length(show_inds),1);
for jj = 1:length(show_inds)
    total_V_explained_integral(jj) = sum(data_params.snp_h_liab(trait_inds) .* ...
        data_params.Correction(trait_inds,jj));
end
total_V_explained_str = ... % show var. explained in %
    num2str_cell(num2cell(total_V_explained_integral*100),3,[],1);
total_h_explained_str = ... % show var. explained in %
    num2str_cell(num2cell(total_V_explained_integral*100 ./ data.h_familial(trait_inds(1))), 3,[],1);
for jj = 1:length(show_inds)
    cur_legend_vec{jj} = [legend_vec{show_inds(jj)} ...
        ', V-expl.=' total_V_explained_str{jj} ', h^2-expl.=' total_h_explained_str{jj}];
end

legend(cur_legend_vec);
switch data.trait_type{trait_inds(1)}
    case {'QTL', 'Quantitative'}
        xlim([0 0.3]); % temp. effect size range ..
    case {'Binary', 'binary'}
        xlim([1 max(data.OR(trait_inds,2))+1]);
end
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'var_explained_by_effect_size_smoothed_density'), 'jpg');
%    my_saveas(gcf, fullfile(figs_dir, trait_is, ...
%        [trait_is '_var_density']), format_fig_vec);

figure; hold on; % plot cummulative var explained as dependency on beta
%                        plot(beta_grid, cumsum(all_beta_estimate));
for jj=1:length(show_inds) % num_corrections
    %        if(show_inds(jj) <= num_corrections)
    plot(beta_grid(end:-1:1), (beta_grid(2)-beta_grid(1)) .* ...
        cumsum(all_beta_estimate_corrected(show_inds(jj),end:-1:1)), ...
        color_vec(jj), 'linewidth', line_width_vec(jj));
end
title([' Cum. Var. explained as function of effect size, ' str2title(disease_names{i}) ]);
xlabel(effect_str); ylabel('Var. Explained'); % should be '|\beta|'
legend(cur_legend_vec);
%                        legend(['observed ' total_V_explained_str{1}], ...
%                            ['corrected ' total_V_explained_str{2}]);
switch data.trait_type{trait_inds(1)}
    case {'QTL', 'Quantitative'}
        xlim([0 0.3]);
    case {'binary', 'Binary'}
        xlim([1 max(data.OR(trait_inds,2))+1]);
end
my_saveas(gcf, fullfile(root_disease_dir, disease_names{i}, ...
    'var_explained_by_effect_size_smoothed_cumulative'), 'jpg');
%    my_saveas(gcf, fullfile(figs_dir, trait_is, ...
%        [trait_is '_var_cumulative']), format_fig_vec);


