% Plot all figures we need for the paper


figs_dir = '../../common_disease_model/docs/pnas/power_paper/figs/automatic';
tables_dir = '../../common_disease_model/docs/pnas/power_paper/docs/tables/automatic';
demographic_sims_dir = '../../common_disease_model/data/schaffner_simulations/EuropeFixed/files_europe_fixed/files_var_eric/mat'; % new file (corrected Europe) !!
equilibrium_parameters_output_file =  '../../common_disease_model/figs/EyreWalker/new_eric/equilibrium/two_class_equilibrium_parameters';
power_figs_dir = fullfile(figs_dir, 'power');



two_class_output_file = fullfile(figs_dir, 'figs_data_two_class_model.mat');
figure_inds = [1 2 3 4 5]; % which figures do we want to display
supp_figure_inds = []; % [1 2 3 4 5]; % which supp. figs. we want to display
table_inds = [1 2 3]; % which tables to save

figs_for_paper_flag = 0; % 1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Parameters:
N = 10000; % effective population size
p_val_cutoff_vec = [0.05 0.05 / 20000];
p_val_cutoff = 0.05 / 20000; % account for # of genes in the genoms
alpha_vec = 1/3; % 0.999999999; % proportion of functional 'null' rare mutations
full_enrichment_alpha_vec = [alpha_vec 1]; % represent how many of alleles are null
beta_vec = linspace(0.001, 10.001, 201); % logspace(-6,1,100); % effect size
mu = 0.01; % prevalence
n_cases_vec = 1000:1000:100000; n_controls_vec = n_cases_vec; n_vec = n_cases_vec + n_controls_vec;
[~, i_beta] = min(abs(beta_vec  - 3)); % 61; % choose specific beta
[~, i_n] = min(abs(n_vec - 50000)); %  25; % choose specific n
prevalence = 0.05; % prevalence (assume 5%)

rare_cumulative_per_gene = 1; % cumulative allele frequency of rare alleles per gene
f_rare = 0.01; % frequency below which an allele is considered 'rare'
f_rare_vec = logspace(-6,0,1000); % different possible maximal values
f_rare_vec(end) = 0.99999999999999;
s_null_vec = [0 logspace(-6,-1,500)]; % try different selection coefficients
show_s_null = [0 logspace(-5, -1, 9)]; % [0 logspace(-4, -2, 5)]; %     [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01];

num_null_s = length(show_s_null);
for i=1:length(show_s_null)
    [~, show_s_null_ind(i)] = min(abs(s_null_vec - show_s_null(i)));
end
s_null = 0.01; % selection coefficient on null (functional) alleles - pretty strong


L = 2000; % number of loci (gene length)
mutation_rate = 2*10^(-8); % mutation rate per-nucleotide per-generation
trait_type = 'quantitative'; % 'disease'; % simulate either disease or quqantitative traits
rare_cumulative_per_gene = 1; % cumulative allele frequency of rare alleles per gene (theta)
theta = 4*N*mutation_rate; % Effective mutation rate
x_vec = (1:2*N-1)./(2*N); % vector of allele frequencies

% Parameters of different types of mutations
mutation_str = {'missense', 'missense-strong', 'missense-neutral', 'stop', 'frameshift', 'synonymous', 'total-substitution'}; %  'synonomous'};
mutation_type_str = {'(mixed)', '(null)', '(neutral)', '(null)', '(null)', '(neutral)', '(mixed)'};
type_str = {'strong-missense', 'neutral-missense'};
mutation_fraction_birth_vec = [2/3-1/20 1/20 1/10 1/3 1]; % The different mutation rates
strong_prop =1/3;
max_f_rare_vec = [0.001 0.01]; % thresholds for being considered rare

%n = 2000; % number of individuals
% prevalence = 0.1; % disease prevalence for disease traits


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demographic_models_struct = [];
demographic_models_struct.file_names = GetFileNames(fullfile(demographic_sims_dir, 'cumul*.mat'), 1);  % NEW! add before/after bottle-neck


for i=1:length(demographic_models_struct.file_names)
    demographic_models_struct.model_str{i} = ...
        strdiff(strdiff(remove_suffix_from_file_name(remove_dir_from_file_name(demographic_models_struct.file_names{i})), ...
        'cumul_'), '_all');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist(two_class_output_file, 'file')) % first look if data already computed
    load(two_class_output_file);
end
if(~exist('w_x_null_mat', 'var'))  % no need to run again and again
    w_x_null_mat = cell(2,1); w_x_harmless = cell(2,1); w_all = cell(2,1);
    c_cumulative = cell(2,1); frac_null_by_freq_cumulative = cell(2,1);
    for j=1:1 % loop on ONE possible mixture
        compute_alpha_j = j
        [two_class_stat_struct ...
            w_x_null_mat{j} w_x_harmless{j} w_all{j} c_cumulative{j} frac_null_by_freq_cumulative{j}] = ...
            compute_two_class_model_parameters(s_null_vec, ...
            f_rare_vec, full_enrichment_alpha_vec(j), rare_cumulative_per_gene, N, two_class_output_file);
    end
    save(two_class_output_file, 'two_class_stat_struct', ...
        'w_x_null_mat', 'w_x_harmless', 'w_all', 'c_cumulative', 'frac_null_by_freq_cumulative');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=figure_inds % create main figures
    switch i        
        case 1 % Figure 1: distribution of rare alleles.
            plot_bayes_factor = 0; % display plots independently of mixture coefficients
            
            title_str = ['Detect. power rare. \alpha=' num2str(alpha_vec*100,3) '%' ...
                ' functional, f^*=' num2str(f_rare*100,2) '%, c=' ...
                num2str(rare_cumulative_per_gene*100,3) '% cum. freq., p-val cutoff = ' num2str(p_val_cutoff, 3)];
            
            plot_two_class_equilibrium_statistics(two_class_stat_struct, w_x_null_mat, w_x_harmless, w_all, ...
                f_rare_vec, frac_null_by_freq_cumulative, c_cumulative, ...
                N, show_s_null, show_s_null_ind, s_null_vec, rare_cumulative_per_gene, alpha_vec, prevalence, ...
                title_str, figs_dir, figs_dir, plot_bayes_factor, figs_for_paper_flag, ...
                demographic_models_struct, equilibrium_parameters_output_file); % Plot properties of distribution - this should also be part of the figures for paper

            
            
            
        case 2 % Figure 2: number of alleles as function of s 
            
            
    end % switch i
end % loop on i

plot_two_class_power_statistics(equilibrium_parameters_output_file, power_figs_dir, power_figs_dir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=supp_figure_inds % create supplamentary figures
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist(two_class_output_file, 'file')) % first look if data already computed
    load(two_class_output_file);
end
for i=table_inds % create tables
    switch i % which table to generate
        case 1 % table 1: #of alleles for different allele frequencies                        
            R = compute_two_class_table(two_class_stat_struct, ...
                mutation_fraction_birth_vec, f_rare_vec, max_f_rare_vec, c_cumulative, ...
                frac_null_by_freq_cumulative, ...
                strong_prop, mutation_str, mutation_type_str, mutation_rate, ...
                mu, theta, L, N, num_null_s, show_s_null, show_s_null_ind);
            
            my_mkdir(tables_dir);
            savecellfile(R, fullfile(tables_dir, 'table1.txt')); % save table to file
    end % switch which table to generate
end % loop on table inds




% for j=1:length(s_null_vec) % loop on different selection coefficients
%     Z_empiric(j) = absorption_time_by_selection(s_null_vec(j), 1, N, 0, f_rare_vec(end), 'freq');
%     Z_theoretic(j) = phi_s_integral(f_rare_vec(end), -4*N*s_null_vec(j), 1) - phi_s_integral(0.000000001, -4*N*s_null_vec(j), 1);
% end
% figure;  semilogx(s_null_vec, Z_empiric, 'linewidth', 4); hold on; semilogx(s_null_vec, Z_theoretic, 'r--', 'linewidth', 2)
%figure; plot(Z_empiric, Z_theoretic)

s = 0.001;  max_f = 0.999;
Z_disease = phi_s_integral(max_f, -4*N*s, 'disease', 1) - ...
    phi_s_integral(0.001, -4*N*s, 'disease', 1)
Z_quant = phi_s_integral(max_f, -4*N*s, 1) - ...
    phi_s_integral(0.001, -4*N*s, 1)

