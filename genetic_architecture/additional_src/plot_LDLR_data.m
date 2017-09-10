% New - plot for LDLR. Can run the code in this function straight from the command-line
function plot_LDLR_data(new_figs_dir)

AssignGeneralConstants;

if(~exist('new_figs_dir', 'var') || isempty(new_figs_dir))
    new_figs_dir = '..\..\common_disease_model\figs\EyreWalker\new_eric';
end
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

n_vec = [5000 500 50];

n_ctr=1;
for n_ctr=1:length(n_vec)  % loop on sample size
    demographic_models_struct.data{1}.rho_vec_finite{n_ctr} = ...
        frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
        10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
        demographic_models_struct, 'europ', n_vec(n_ctr)) % use new function - finite sample
    
    
    demographic_models_struct.data{1}.rho_vec_finite_sampling{n_ctr} = ...
        frac_null_conditional_on_freq_less_f( ...
        10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
        demographic_models_struct, 'europ', n_vec(n_ctr), [], 'sampling', 20000*5000/n_vec(n_ctr)) % use new function - finite sample with sampling!!!
end

% figure; semilogx(x_vec, demographic_models_struct.data{1}.f_null_vec(1, :), 'r'); hold on; % null
% semilogx(x_vec, demographic_models_struct.data{1}.f_null_vec(2, :), 'b'); hold on; % neutral

figure; semilogx(x_vec, demographic_models_struct.data{1}.rho_vec, 'k', 'linewidth', 2); hold on; % rho vec population



% semilogx(x_vec, demographic_models_struct.data{1}.rho_vec2, 'b--');

for n_ctr = 1:length(n_vec)
    first_ind = find(demographic_models_struct.data{1}.rho_vec_finite{n_ctr} < alpha, 1);
    semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite{n_ctr}(first_ind:end), ...
        color_vec(n_ctr), 'linewidth', 2);
    semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite_sampling{n_ctr}(first_ind:end), ...
        [color_vec(n_ctr) '--'], 'linewidth', 2);
end


%semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite2(first_ind:end), 'm', 'linewidth', 2);
%semilogx(x_vec(first_ind:end), demographic_models_struct.data{1}.rho_vec_finite_sampling(first_ind:end), 'g', 'linewidth', 2);

%legend({'population', 'finite-sample'}); % , 'finite-analytic2', 'finite-sampling'});
%legend({'population', 'finite-analytic', 'finite-analytic2', 'finite-sampling'});

legend_vec =  cellstr([repmat('n=', length(n_vec), 1) num2str(n_vec')]);
legend_vec = ['pop-freq', mat2vec([legend_vec legend_vec]')']; 
for i=1:length(n_vec)
    legend_vec{2*i+1} = [legend_vec{2*i+1} ' (sim.)']; 
end
% legend_vec = {'pop-freq', 
h_leg = legend(legend_vec, 1);
set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
xlim([10^(-6), 1]);
xlabel('Derived allele frequency, f'); ylabel('Proportion of nulls, \rho_s(f)');
add_faint_grid(0.5);
my_saveas(gcf, fullfile(new_figs_dir, 'rho_finite_sample'), {'epsc', 'pdf', 'jpg'});

n_vec = round(logspace(1,5, 100)); 
%n_vec = [10 50 100 200 500 1000 2000 5000 10000];
for n_ctr=1:length(n_vec)  % loop on sample size
    run_n = n_vec(n_ctr)
    rho_by_sample_size(n_ctr) = ...
        frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
        10^(-1.7), [], alpha, max(10^(-4), 1/(2*n_vec(n_ctr))), ...
        demographic_models_struct, 'europ', n_vec(n_ctr)); % use new function - finite sample
end

rho_population = ...
    frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
    10^(-1.7), [], alpha, 10^(-4), ...
    demographic_models_struct, 'europ'); % use new function - finite sample

figure; semilogx(n_vec, rho_by_sample_size, 'linewidth', 2);     % New: plot different sample sizes
hold on; 
plot(n_vec, repmat(rho_population(220), length(n_vec), 1), 'r--', 'linewidth', 2); 

h_leg = legend({'finite-sample', 'population'}, 4); 
set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
xlabel('Number of cases'); ylabel('Proportion of nulls, \rho_s(f)');
add_faint_grid(0.5); 
my_saveas(gcf, fullfile(new_figs_dir, 'rho_finite_sample_singletons'), {'epsc', 'pdf', 'jpg'});


% Get values for two thresholds we seek 
T_vec = [0.0001 0.01]; % all
for j=1:length(T_vec)
    [~, T_inds(j)] = min(abs(x_vec - T_vec(j)));
end

% % % rho_ideal = demographic_models_struct.data{1}.rho_vec(T_inds)
% % % rho_finite = demographic_models_struct.data{1}.rho_vec_finite(T_inds)'
% % % rho_polyphen_ideal = demographic_models_struct.data{1}.rho_vec_polyphen(T_inds)
% % % rho_polyphen_finite = demographic_models_struct.data{1}.rho_vec_polyphen_finite(T_inds)
% % % 
% % % 
% % % lambda_ideal = lambda .* demographic_models_struct.data{1}.rho_vec(T_inds)
% % % lambda_finite = lambda .* demographic_models_struct.data{1}.rho_vec_finite(T_inds)'
% % % lambda_polyphen_ideal = lambda .* demographic_models_struct.data{1}.rho_vec_polyphen(T_inds)
% % % lambda_polyphen_finite = lambda .* demographic_models_struct.data{1}.rho_vec_polyphen_finite(T_inds)



% % % Additional tries, including polyphen!! 
% % % demographic_models_struct.data{1}.rho_vec_finite2 = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
% % %     10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
% % %     demographic_models_struct, 'europ', 5000) % use new function - finite sample
% % % demographic_models_struct.data{1}.rho_vec_finite_sampling500 = frac_null_conditional_on_freq_less_f( ...
% % %     10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
% % %     demographic_models_struct, 'europ', 500, [], 'sampling') % use new function - finite sample with sampling!!!
% % % demographic_models_struct.data{1}.rho_vec_finite_sampling50 = frac_null_conditional_on_freq_less_f( ...
% % %     10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
% % %     demographic_models_struct, 'europ', 50, [], 'sampling') % use new function - finite sample with sampling!!!
% % % demographic_models_struct.data{1}.rho_vec_polyphen_finite = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
% % %     10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
% % %     demographic_models_struct, 'europ', 5000, [0.2 0.2]) % use new function - finite sample
% % % demographic_models_struct.data{1}.rho_vec_polyphen2 = frac_null_conditional_on_freq_less_f(  ... % NOT working yet !!
% % %     10^(-1.7), [], alpha, demographic_models_struct.data{1}.x_vec, ...
% % %     demographic_models_struct, 'europ', [], [0.2 0.2]) % use new function - finite sample
