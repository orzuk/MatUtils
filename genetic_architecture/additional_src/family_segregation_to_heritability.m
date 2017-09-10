% New - plot for LDLR. Can run the code in this function straight from the command-line
function plot_LDLR_data()


% NEED TO DEBUG!!!

demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/EuropeFixed/files_europe_fixed/files_var_eric'; % New Version !! Europe fixed
figs_dir = '../../common_disease_model/figs/EyreWalker/new_eric/all_models/empirical_distribution';
allele_freq_file = fullfile(demographic_sims_dir, 'fbin.txt');

XXX = load('../../common_disease_model/data/schaffner_simulations/EuropeFixed/fine_s_var/out_var/mat/europ.1.7.mat');
XXX_NEUTRAL = load('../../common_disease_model/data/schaffner_simulations/EuropeFixed/fine_s_var/out_var/mat/europ.0.0.mat');

XXX.s_cumulative = cumsum(XXX.s_distribution, 2);
XXX_NEUTRAL.s_cumulative = cumsum(XXX_NEUTRAL.s_distribution, 2);


lambda = 17.1; % effect size
alpha = 0.25; % fraction of nulls

num_bins = 500; smooth_width = 5;
x_vec = ReadDataFile(allele_freq_file);

demographic_models_struct = [];
demographic_models_struct.model_str{1} = 'europ';
demographic_models_struct.data{1}.x_vec = x_vec.derived_freq; % Note: allele freq. is in logarithmic coordinates.
demographic_models_struct.data{1}.p_vec(2,:) = mean(XXX_NEUTRAL.s_cumulative);
demographic_models_struct.data{1}.p_vec(1,:) = mean(XXX.s_cumulative);
demographic_models_struct.data{1}.num_alleles_per_chrom = demographic_models_struct.data{1}.p_vec(:,end);
demographic_models_struct.data{1}.rho_vec = demographic_models_struct.data{1}.p_vec(1,:) .* alpha ./ ...
    (demographic_models_struct.data{1}.p_vec(1,:) .* alpha + demographic_models_struct.data{1}.p_vec(2,:) .* (1-alpha));
demographic_models_struct.data{1}.f_null_vec = bsxfun(@rdivide, demographic_models_struct.data{1}.p_vec(end:-1:1,:), ...
    demographic_models_struct.data{1}.p_vec(end:-1:1,end));  % make sure this is normalized. This should go UP!
x_vec = demographic_models_struct.data{1}.x_vec;

demographic_models_struct.data{1}.s_vec = [10^(-1.7) 0];

demographic_models_struct.data{1}.rho_vec = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ') % use new function - finite sample
demographic_models_struct.data{1}.rho_vec_polyphen = demographic_models_struct.data{1}.p_vec(1,:) .* alpha .* 0.8 ./ ...
    (demographic_models_struct.data{1}.p_vec(1,:) .* alpha .* 0.8 + demographic_models_struct.data{1}.p_vec(2,:) .* (1-alpha) .* 0.2);


demographic_models_struct.data{1}.rho_vec_finite = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ', 5000) % use new function - finite sample

demographic_models_struct.data{1}.rho_vec_finite2 = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ', 5000) % use new function - finite sample

demographic_models_struct.data{1}.rho_vec_finite_sampling = frac_null_conditional_on_freq_less_f( ...
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ', 15000, [], 'sampling') % use new function - finite sample with sampling!!!

demographic_models_struct.data{1}.rho_vec_finite_sampling500 = frac_null_conditional_on_freq_less_f( ...
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ', 500, [], 'sampling') % use new function - finite sample with sampling!!!

demographic_models_struct.data{1}.rho_vec_finite_sampling50 = frac_null_conditional_on_freq_less_f( ...
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ', 50, [], 'sampling') % use new function - finite sample with sampling!!!

demographic_models_struct.data{1}.rho_vec_polyphen_finite = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ', 5000, [0.2 0.2]) % use new function - finite sample

demographic_models_struct.data{1}.rho_vec_polyphen2 = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
    10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
    demographic_models_struct, 'europ', [], [0.2 0.2]) % use new function - finite sample

figure; semilogx(x_vec, demographic_models_struct.data{1}.f_null_vec(1, :), 'r'); hold on; % null
semilogx(x_vec, demographic_models_struct.data{1}.f_null_vec(2, :), 'b'); hold on; % neutral

figure; semilogx(x_vec, demographic_models_struct.data{1}.rho_vec, 'b', 'linewidth', 2); hold on; % rho vec
% semilogx(x_vec, demographic_models_struct.data{1}.rho_vec2, 'b--');
first_ind = find(demographic_models_struct.data{1}.rho_vec_finite < alpha, 1);
semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite_sampling(first_ind:end), 'r', 'linewidth', 2);
semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite_sampling500(first_ind:end), 'g', 'linewidth', 2);
semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite_sampling50(first_ind:end), 'm', 'linewidth', 2);


%semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite2(first_ind:end), 'm', 'linewidth', 2);
%semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite_sampling(first_ind:end), 'g', 'linewidth', 2);

%legend({'population', 'finite-sample'}); % , 'finite-analytic2', 'finite-sampling'});
%legend({'population', 'finite-analytic', 'finite-analytic2', 'finite-sampling'});
h_leg = legend({'pop-freq.', 'finite-sample'}, 1);
set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
xlim([10^(-6), 1]);
xlabel('Derived allele frequency, f'); ylabel('Proportion of nulls, \rho_s(f)');
add_faint_grid(0.5);
my_saveas(gcf, fullfile(new_figs_dir, 'rho_finite_sample'), {'epsc', 'pdf', 'jpg'});

figure; semilogx(x_vec, demographic_models_struct.data{1}.rho_vec ./ ...
    demographic_models_struct.data{1}.rho_vec_finite', 'b', 'linewidth', 2); hold on; % rho vec



T_vec = [0.0001 0.01]; % all
for j=1:length(T_vec)
    [~, T_inds(j)] = min(abs(x_vec - T_vec(j)));
end

rho_ideal = demographic_models_struct.data{1}.rho_vec(T_inds)
rho_finite = demographic_models_struct.data{1}.rho_vec_finite(T_inds)'
rho_polyphen_ideal = demographic_models_struct.data{1}.rho_vec_polyphen(T_inds)
rho_polyphen_finite = demographic_models_struct.data{1}.rho_vec_polyphen_finite(T_inds)


lambda_ideal = lambda .* demographic_models_struct.data{1}.rho_vec(T_inds)
lambda_finite = lambda .* demographic_models_struct.data{1}.rho_vec_finite(T_inds)'
lambda_polyphen_ideal = lambda .* demographic_models_struct.data{1}.rho_vec_polyphen(T_inds)
lambda_polyphen_finite = lambda .* demographic_models_struct.data{1}.rho_vec_polyphen_finite(T_inds)

