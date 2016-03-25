% Script for runnin sfs_code forward simulator for Fisher-Wright Model

AssignGeneralConstants;
orange = [1 0.6 0.1];
population_color_vec = {'k', orange, 'r', 'g', 'b', 'c'};
% eric_color_vec = 'kbgyr'; % Conenstion for selection coefficients (we don't have orange. Use yellow)
eric_color_vec = {'k', 'b', 'g', orange, 'r'}; % Conenstion for selection coefficients (we don't have orange. Use yellow)

my_symbol_vec = {'--', '-', ':', '-.'}; % flip ordering (set integer powers as solid lines)
selection_color_vec = {'k', 'b', 'g', orange, 'r'}; % replace yellow with orange
model_str = 'equilibrium';
two_side_flag = 0; % 0 - use derived allele frequency, 1 - use minor allele frequency
sim_flag = 'mine'   % sfs , schaffner


if(machine == UNIX) % set sfs_code path
    % sfs_code_dir = 'c/Users/user/Downloads/sfscode';
else % PC
    sfs_code_dir = 'C:/Users/user/Downloads/sfscode';
end

s_vec =  [10.^(-1:-0.5:-5) 0]; % , s=10^-1.5, s=10^-2, ... s=10^-5, s=0);
num_s = length(s_vec);
s_vec_str = cellstr([repmat('10^{', num_s, 1) num2str(log10(vec2column(s_vec)), 2) repmat('}', num_s, 1)]);
for i=find(s_vec == 0)
    s_vec_str{i} = '0';
end
s_vec_str = strrep_cell(s_vec_str, ' ', '');
for i=1:length(s_vec_str)
    s_vec_str{i} = ['s=' s_vec_str{i}];
    
    sfs_output_gazava{i} = fullfile(sfs_code_dir, 'scripts', 'Tennenssen.out'); % update here according to s 
end

sfs_sims_dir = '../../common_disease_model/data/schaffner_simulations/EuropeFixed/files_europe_fixed/files_var_eric'; % New Version !! Europe fixed
sfs_figs_dir = 'C:\Users\user\Dropbox\rare_alleles_paper\JamesZou\figs';

% set demoraphic model: convention:
iters = 1000; % number of iterations
num_generations = 100; % [];
N_vec = [10000 50000]; % [];  % effective population size
%N_vec = ceil(10.^(2.5:0.5:6));
num_N = length(N_vec);

mu = 10^(-5); % mutation rate
expansion_factor = 1; % no expansion
cum_flag = 'weighted';


% Create models for dempgraphy
% Set also other demographic_models
N_vec_gazava = create_demographic_model('gazava');


% Simulations ! takes a long time
for i=1:length(s_vec) % loop on selection
    run_s = s_vec(i)
    num_generations_gazava =  length(N_vec_gazava)-1;
    
    switch lower(sim_flag)
        case 'mine'  % run simulations using own Matlab program
            [freq_struct_gazava{i}, absorption_struct_gazava{i}, simulation_struct_gazava{i}, ...
                N_vec_out_gazava{i}, simulation_time_gazava{i}] = ... % New: separate output to different structures
                FisherWrightSimulation(N_vec_gazava, mu, -abs(s_vec(i)), num_generations_gazava, expansion_factor, ...
                'equilibrium', iters, 'simulation')
            x_vec_gazava{i} = freq_struct_gazava{i}.x_vec{num_generations_gazava}./max(freq_struct_gazava{i}.x_vec{num_generations_gazava});
            y_vec_gazava{i} = cumsum(freq_struct_gazava{i}.p_vec{num_generations_gazava} .* freq_struct_gazava{i}.x_vec{num_generations_gazava});
            y_vec_gazava{i} = y_vec_gazava{i} ./ y_vec_gazava{i}(end);  % normalize
            y_vec_not_weighted_gazava{i} = cumsum(freq_struct_gazava{i}.p_vec{num_generations_gazava}) ./ ...
                sum(freq_struct_gazava{i}.p_vec{num_generations_gazava});
        case 'sfs' % read data from sfs-code simulations
            freq_struct_gazava{i} = parse_sfs_code_output(sfs_output_gazava{i}); % read output
            x_vec_gazava{i} = freq_struct_gazava{i}.freq_vec;
            y_vec_gazava = cumsum(x_vec_gazava{i});
            y_vec_gazava{i} = y_vec_gazava{i} ./ y_vec_gazava{i}(end);  % normalize
            y_vec_not_weighted_gazava{i}
            
        case 'schaffner' % run simulations using Schaffner's program 
            
    end % switch sim_flag    
end % loop on selection coefficient s 

N_vec = ceil(10.^(2.5:0.5:6));
for i=1:length(s_vec) % loop on selection
    for i_N = 4 % 1:length(N_vec)  % loop on population size
        run_N = N_vec(i_N)
        run_s2 = s_vec(i)
        [freq_struct{i,i_N}, absorption_struct{i,i_N}, simulation_struct{i,i_N}, ...
            N_vec_out{i,i_N}, simulation_time{i,i_N}] = ... % New: separate output to different structures
            FisherWrightSimulation(N_vec(i_N), mu, -abs(s_vec(i)), num_generations, expansion_factor, ...
            'equilibrium', iters, 'simulation')
        x_vec{i,i_N} = freq_struct{i,i_N}.x_vec{num_generations}./max(freq_struct{i,i_N}.x_vec{num_generations});
        y_vec{i,i_N} = cumsum(freq_struct{i,i_N}.p_vec{num_generations} .* freq_struct{i,i_N}.x_vec{num_generations});
        y_vec{i,i_N} = y_vec{i,i_N} ./ y_vec{i,i_N}(end);  % normalize
        y_vec_not_weighted{i,i_N} = cumsum(freq_struct{i,i_N}.p_vec{num_generations}) ./ ...
            sum(freq_struct{i,i_N}.p_vec{num_generations});
    end % loop on N
end % loop on S (selection)

%%%%%%%%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = 0.0001;
i_N = 4; % set N=10,000
final_x_vec = (1:2*N_vec(i_N)-1) ./ (2*N_vec(i_N)); % set new coordinates
for i=1:length(s_vec) % loop on selection
    tmp_p_vec_equlibrium_analytic{i} = exp( allele_freq_spectrum([0 final_x_vec 1], -s_vec(i), N_vec(i_N), two_side_flag, 'log', 1) ); % compute analytic approxiamtion (valid only for constant population size)
    tmp_p_vec_equlibrium_analytic2{i} = phi_s_integral([0 final_x_vec 1], -s_vec(i)*4*N_vec(i_N), -1 ); % compute analytic approxiamtion (valid only for constant population size)
    
    tmp_p_vec_equlibrium_analytic{i} = normalize_hist(final_x_vec, tmp_p_vec_equlibrium_analytic{i}(2:end-1)); % normalized
    
    mean_weighted_freq(i) = absorption_time_by_selection(s_vec(i), theta, N_vec(1), 1/(2*N_vec(i_N)), 1- 1/(2*N_vec(i_N)), 2) / ...
        absorption_time_by_selection(s_vec(i), theta, N_vec(i_N), 1/(2*N_vec(1)), 1- 1/(2*N_vec(i_N)), 1);
    
    tmp_median_equilibrium_analytic(i) = final_x_vec(find(cumsum(tmp_p_vec_equlibrium_analytic{i}) / ...
        sum(tmp_p_vec_equlibrium_analytic{i}) >= 0.5, 1));
end


i_N = 4; % set N=10,000
final_x_vec = (1:2*N_vec(i_N)-1) ./ (2*N_vec(i_N)); % set new coordinates
for plot_simulations = 2 % 0:2
    figure; % hold on;
    for i=1:length(s_vec) % loop on selection
        tmp_color_ind = mod_max(1+floor(mod(10-i, 10)/2), 5);
        
        % compute and plot simulation !!!
        switch plot_simulations
            case 2 % simulation expansion gazave model
                model_str = 'gazava';
                sim_str = '_expansion_gazava';
                semilogx(x_vec_gazava{i}, y_vec_gazava{i}, ...
                    'color', eric_color_vec{tmp_color_ind}, ...
                    'linestyle', my_symbol_vec{0+mod_max(i-1,2)}, 'linewidth', 2); hold on;  % plot analytic
                
            case 1  % simulation equilibrium
                model_str = 'equilibrium';
                sim_str = '_and_simulations';
                semilogx(x_vec{i,i_N}, y_vec{i,i_N}, ...
                    'color', eric_color_vec{tmp_color_ind}, ...
                    'linestyle', my_symbol_vec{0+mod_max(i-1,2)}, 'linewidth', 2); hold on;  % plot analytic
                
            case 0
                model_str = 'equilibrium';
                semilogx(final_x_vec, cumsum(tmp_p_vec_equlibrium_analytic{i})/sum(tmp_p_vec_equlibrium_analytic{i}), ...
                    'color', eric_color_vec{tmp_color_ind}, ...
                    'linestyle', my_symbol_vec{mod_max(i-1,2)}, 'linewidth', 2); hold on;  % plot analytic
                sim_str = '_analytic';
        end
    end
    xlabel(str2title('f')); ylabel('Cumulative Freq.');
    ylim([0, 1.01]);
    switch plot_simulations
        case 2
            xlim([10^(-7) 1]);
        otherwise
            xlim([10^(-4.32) 1]);
    end % switch simulation type
    h_leg = legend(s_vec_str, 'location', 'northwest'); % ([end:-1:2 1])
    title(['Variation in cumulative allele frequency. ' model_str], 'fontsize', 14, 'fontweight', 'bold'); %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
    add_faint_grid(0.5);
    set(h_leg, 'edgecolor',[0.99 0.99 0.99]); % 'Ycolor',[0.8 0.8 0.8]);
    my_saveas(gcf, fullfile(sfs_figs_dir, ['cumulative_freq_' model_str  sim_str]), {'epsc', 'pdf'}); % save file
end % loop on plot simulations


% Just check that resulting medians make sense
%figure; semilogx(s_vec, mean_weighted_freq); xlabel('s'); ylabel('mean-weighted-freq'); title('Mean Weighted Frequency');
figure; semilogx(s_vec, tmp_median_equilibrium_analytic); xlabel('s'); ylabel('meadin-weighted-freq'); title('Median Weighted Frequency');


% Compute mean and median for each s
s_null_vec = [10.^(-1:-0.01:-5) 0]; % , s=10^-1.5, s=10^-2, ... s=10^-5, s=0);
N_vec = ceil(10.^(2.5:0.5:6)); num_N = length(N_vec);

num_s_null = length(s_null_vec);
mean_x = zeros(num_N,num_s_null);
median_x = zeros(num_N,num_s_null); median_mixture_x = zeros(num_N,num_s_null);
mean_het_x = zeros(num_N,num_s_null); total_het_x = zeros(num_N,num_s_null);
integral_phi_x = zeros(num_N,num_s_null);
normalization_factor_x = zeros(num_N,num_s_null);

for j_N = 1:length(N_vec) % loop on N !
    N=N_vec(j_N); % 10000;
    S_vec = s_null_vec .* N .* 4;
    for i=1:length(s_null_vec)
        S = s_null_vec(i) * N * 4;
        integral_phi_x(i) = phi_s_integral(0.999999999,-S, 0) - phi_s_integral(0.000000001,-S, 0);
        mean_x(j_N,i) = phi_s_integral(0.999999999,-S, 2) - phi_s_integral(0.000000001,-S, 2);
        normalization_factor_x(j_N,i) = phi_s_integral(0.999999999,-S, 1) - phi_s_integral(0.000000001,-S, 1);
        mean_x(j_N,i) = mean_x(j_N,i) ./ normalization_factor_x(j_N,i); % normalize to get mean of polymorphic alleles
        median_x(j_N,i) = fzero(@(x) phi_s_integral(x,-S, 1) - ...
            0.5 .* phi_s_integral(0.999999999,-S, 1) - 0.5*phi_s_integral(0.000000001,-S, 1), mean_x(j_N,i)); % solution of phi(x)-phi(0) = (1/2)*( phi(1)-phi(0))
        mean_het_x(j_N,i) = phi_s_integral(0.999999999,-S, -2) - phi_s_integral(0.000000001,-S, -2);
        mean_het_x(j_N,i) = mean_het_x(j_N,i) ./  normalization_factor_x(j_N,i);
        total_het_x(j_N,i) = phi_s_integral(0.999999999,-S, 'var') - phi_s_integral(0.000000001,-S, 'var');
        i_is = i
    end
end % loop on N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% FIGURE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop on N !
for median_flag = 0:1
    if(median_flag)
        median_str = 'Median';
    else
        median_str = 'Mean';
    end
    figure;
    for i_N = 1:length(N_vec)
        if(median_flag)
            plot_y_vec = median_x(i_N,:);
        else
            plot_y_vec = mean_x(i_N,:);
        end
        %    loglog(s_null_vec, mean_x, 'r', 'linewidth', 2); hold on;
        loglog(s_null_vec, plot_y_vec, ...
            'color', eric_color_vec{ceil(i_N/2)}, ...
            'linestyle', my_symbol_vec{mod_max(i_N,2)}, ...
            'linewidth', 2); hold on; %'b'
    end
    xlabel('s'); ylabel('f'); title([median_str ' frequencies as function of s and N']);
    
    add_faint_grid(0.5);
    N_vec_str = cellstr([repmat('10^{', num_N, 1) num2str(log10(vec2column(N_vec)), 2) repmat('}', num_N, 1)]);
    for i=find(N_vec == 0)
        N_vec_str{i} = '0';
    end
    N_vec_str = strrep_cell(N_vec_str, ' ', '');
    for i=1:length(N_vec_str)
        N_vec_str{i} = ['N=' N_vec_str{i}];
    end
    %    h_leg = legend({'Mean', 'Median'});
    h_leg = legend(N_vec_str, 'location', 'southwest');
    set(h_leg,'edgecolor',[0.99 0.99 0.99]); % ,'Ycolor',[0.8 0.8 0.8]);
    my_saveas(gcf, fullfile(sfs_figs_dir, [median_str '_freq_' model_str '_analytic']), {'epsc', 'pdf'}); % save file
end % loop on median/mean





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Figure 2.5 %%%%%%%%%%%%%%%%%%%%%%%%

% Compute and plot fixation probability
for j_N = 1:length(N_vec) % loop on N !
    N=N_vec(j_N); % 10000;
    S_vec = s_null_vec .* N .* 4;
    for i=1:length(s_null_vec)
        P_fix(i,j_N) = fixation_prob_by_selection(-s_null_vec(i), N, 1/(2*N));
    end
end
figure;
for i_N = 1:length(N_vec)
    loglog(s_null_vec, P_fix(:,i_N), ...
        'color', eric_color_vec{ceil(i_N/2)}, ...
        'linestyle', my_symbol_vec{mod_max(i_N,2)}, ...
        'linewidth', 2); hold on; %'b'
end
xlabel('s'); ylabel('Pr(fix)'); title('Fixation probability as function of s and N');
ylim([10^(-30), 10^(-2)]);
add_faint_grid(0.5);
N_vec_str = cellstr([repmat('10^{', num_N, 1) num2str(log10(vec2column(N_vec)), 2) repmat('}', num_N, 1)]);
for i=find(N_vec == 0)
    N_vec_str{i} = '0';
end
N_vec_str = strrep_cell(N_vec_str, ' ', '');
for i=1:length(N_vec_str)
    N_vec_str{i} = ['N=' N_vec_str{i}];
end

%    h_leg = legend({'Mean', 'Median'});
h_leg = legend(N_vec_str, 'location', 'northeast');
set(h_leg,'edgecolor',[0.99 0.99 0.99]); % ,'Ycolor',[0.8 0.8 0.8]);

my_saveas(gcf, fullfile(sfs_figs_dir, ['fixation_prob_' model_str '_analytic']), {'epsc', 'pdf'}); % save file



% Plot:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% FIGURE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; % hold on;
for i=1:length(s_vec) % loop on selection
    tmp_color_ind = mod_max(1+floor(mod(i-1, 10)/2), 5);
    semilogx(x_vec{i}, y_vec2{i}, 'color', eric_color_vec{tmp_color_ind}, 'linestyle', ...
        my_symbol_vec{mod_max(i,2)}, 'linewidth', 2); hold on; % what is y_vec2? 
end
xlabel(str2title('f')); ylabel('Cumulative Freq.');
output_plot_file_name = '_num_alleles_dist';
h_leg = legend(s_vec_str([end:-1:2 1]));
title(['Variation in cumulative allele frequency. ' model_str], 'fontsize', 14, 'fontweight', 'bold'); %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
add_faint_grid(0.5);
set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
my_saveas(gcf, fullfile(sfs_figs_dir, ['sfs_' model_str]), {'epsc', 'pdf'}); % save file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% FIGURE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for i=1:length(s_vec) % loop on selection
    tmp_color_ind = mod_max(1+floor(mod(i-1, 10)/2), 5);
    z_vec = cumsum((0.*freq_struct{i}.final_x_vec+1) .* freq_struct{i}.p_vec_equlibrium_analytic);
    z_vec = z_vec ./ z_vec(end);
    semilogx(freq_struct{i}.final_x_vec, z_vec, 'color', eric_color_vec{tmp_color_ind}, ...
        'linestyle', my_symbol_vec{mod_max(i,2)}, 'linewidth', 2); hold on;
end
xlabel(str2title('f')); ylabel('Cumulative Freq. (analytic)');
h_leg = legend(s_vec_str([end:-1:2 1]));
title(['Variation in cumulative allele frequency. (analytic) ' model_str], 'fontsize', 14, 'fontweight', 'bold'); %   Mean=' num2str(mean(f_null)) '. St.d.=' num2str(std(f_null))]);
add_faint_grid(0.5);
set(h_leg,'Xcolor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8]);
my_saveas(gcf, fullfile(sfs_figs_dir, ['sfs_' model_str '_analytic']), {'epsc', 'pdf'}); % save file






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Here do simulations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





