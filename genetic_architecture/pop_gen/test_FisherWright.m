% Script for testing calculations with the Wright-Fisher model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decide which tests to run:
run_sim_flag = 1; % compute site frequency distribution for each s using simulations
test_absorption_time = 0; % compare absorption time for simulations vs. analytical formula
test_moments = 0; % NEW: test analytic moments calculation of SFS using Ewens paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set flags for running
debug_figures = 0;
run_two_expansions_flag = 0; % a mode within simulation, running a two-stage expansion model
one_plot_flag = 1; % make one plot for each s
unite_results_flag = 0; % make one plot for all s values together
save_in_mathematica = 0; % save in mathematica format for Eric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set running parameters
AssignGeneralConstants; AssignRVASConstants;
init_str = 'equilibrium'; % 'equilibrium' 'newly_born'; % start at newly born allele or equilibrium
demography_str = 'two-stage-expan'; % 'expansion1'; % choose demographic model
s_vec = -[0 logspace(-5, -1, 9)]; % -[0 0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01]; %  0.05 0.1]; % -0.00000001; % selection coefficient % s_vec = -[0 logspace(-6, -1, 11)]; % take log-space
compute_mode = 'simulation'; % 'simulation'; % 'simulation';  % 'simulation'; % 'numeric'; % 'simulation'; % 'numeric'; % how to advance calculation
N_vec = []; [N_vec{1}, D] = create_demographic_model(demography_str, 50); N = N_vec{1}(1); % 'expansion1');
mu = mu_per_site * (10000 / N); % mutation rate (per nucleotide per generation)
iters = [20000 20000]; % NEW! always take 2000 % Old: relevant only for simulation. Spend more iterations on equilibrium (new alleles take more time per iteration)
%N = 10; % population size
%num_generations_vec = [200 10]; % two-stage model: slow and fast expansion % 2500; % 50;
%expansion_factor_vec = [1.005 1.05]; % 1.1; % 1.02; % two-stage moe: growth in population size
if(~save_in_mathematica) % dave in mathematica format
    num_bins = 2000; % for plotting histograms
    mathematica_str = '';
    mathematica_flag = 0;
else
    num_bins = [0 logspace(-6, 0, 601)]; % use logarithmic bins  (for Eric)
    mathematica_str = '_mathematica';
    mathematica_flag = 1;
end
fisher_wright_output_dir = '../../common_disease_model/figs/RVAS_gene_specific/FisherWright/';
s_ctr = 1; % counter of selection coefficient
het_struct = cell(length(s_vec), 1);

% New: Set plotting parameters
plot_params.figure_type = 1; plot_params.figs_dir = exome_data_figs_dir; plot_params.hist = 1; plot_params.xlim = [10^(-4) 1];
plot_params.cum=1; plot_params.weighted = 1; plot_params.normalize=1; plot_params.font_size=8; % plot cumulative weighted allele frequency distribution


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run FisherWright Simulation to compute SFS
total_time = cputime;
% Choose model
D.s_grid = s_vec; D.iters = 1000; % take twice as many iters ! 
[D.SFS.x_vec, D.SFS.p_vec] = deal(cell(length(s_vec), 2)); D.SFS.compute_mode = {'simulation', 'numeric'};
for s_ctr = 1:length(s_vec) % 0 % [] %  s_vec(1) %   1:end-2) % (2:end) % s_vec(end) % s_vec(3:end) % (end-1) % s_vec % loop on different selection coefficients    
    % close all; 
    s = s_vec(s_ctr); sprintf('run s=%lf\n', s)
    frac_polymorphic = 2*N*mu * absorption_time_by_selection(s, 1, N, 1/(2*N), 0.999999999, 0);
    for expansion_model = 0:0 % run_two_expansions_flag % allow one rate and two-rate exponential expansions
        ctr=1; % counter on init str
        tmp_dir = demography_str; % new: tmp_dir =
        for init_str = {'equilibrium'} % , 'newly_born'} % , 'equilibrium' 'newly_born'} % list two distributions separately
            fisher_wright_output_file = fullfile(fisher_wright_output_dir, tmp_dir, ['fisher_wright_expansion_' init_str{1}]);
            all_s_output_file = fullfile(fisher_wright_output_dir, tmp_dir, ['all_s_' ...
                init_str{1} mathematica_str]);
            if(run_sim_flag)
                compute_mode_ctr=1; [freq_struct, absorption_struct, simulation_struct] = deal(cell(2,1)); simulation_time = zeros(2,1);
                for compute_mode_ctr = 1:length(D.SFS.compute_mode) %  = {'simulation', 'numeric'}
                    % New!! Go to higher level - compute directly x-vec and p-vec
                    [D.SFS.x_vec{s_ctr, compute_mode_ctr}, D.SFS.p_vec{s_ctr, compute_mode_ctr}, L_correction_factor, compute_time, k_vec, n_vec, weights_vec] = ...
                        compute_allele_freq_spectrum_from_demographic_model(D, s, D.SFS.compute_mode{compute_mode_ctr}, [], mu);
                    %                    compute_allele_freq_spectrum_from_demographic_model(
                    
                    %                    [freq_struct{compute_mode_ctr}, absorption_struct{compute_mode_ctr}, simulation_struct{compute_mode_ctr}, ...
                    %                        N_vec{compute_mode_ctr}, simulation_time(compute_mode_ctr)] = ...
                    %                        FisherWrightSimulation([], D, mu, s, init_str{1}, iters(ctr), compute_mode{1}, num_bins); % run simulation/Markov chain (can take long time!)
                    %compute_mode_ctr=compute_mode_ctr+1;
                end
                my_mkdir(fullfile(fisher_wright_output_dir, tmp_dir));
                save([fisher_wright_output_file '.mat'], ...
                    'D', 'N_vec', 's'); %                    'freq_struct', 'absorption_struct', 'simulation_struct', 'N_vec', 'simulation_time');
            else % already simulated. Load from file
                if(one_plot_flag) % no need to load if not plotting anything
                    load([fisher_wright_output_file '.mat']);
                end
            end
            if(one_plot_flag) % one plot
                % %                 switch init_str{1}
                % %                     case 'equilibrium' % This assumes equilibrium is called first !!!
                % %                         for compute_mode_ctr=1:2
                % %                             save_all_old_x_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.x_vec{D.total_generations}; %  .* (2*N_vec(num_generations));
                % %                             save_all_old_p_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.p_vec{D.total_generations};
                % %                             save_all_old_het_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.het_vec{D.total_generations};
                % %                             save_old_total_het_at_each_generation_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.total_het_at_each_generation_vec
                % %                             if(compute_mode_ctr==1) % simulations
                % %                                 save_old_num_simulated_polymorphic_alleles_vec{compute_mode_ctr} = simulation_struct{compute_mode_ctr}.num_simulated_polymorphic_alleles_vec;
                % %                             end
                % %                         end
                % %                     case 'newly_born'
                % %                         for compute_mode_ctr=1:2
                % %                             freq_struct{compute_mode_ctr}.all_old_x_vec = save_all_old_x_vec{compute_mode_ctr};
                % %                             freq_struct{compute_mode_ctr}.all_old_p_vec = save_all_old_p_vec{compute_mode_ctr};
                % %                             freq_struct{compute_mode_ctr}.all_old_het_vec = save_all_old_het_vec{compute_mode_ctr};
                % %                             freq_struct{compute_mode_ctr}.old_total_het_at_each_generation_vec = save_old_total_het_at_each_generation_vec{compute_mode_ctr};
                % %
                % %                             if(compute_mode_ctr==1) % simulations
                % %                                 simulation_struct{compute_mode_ctr}.new_num_simulated_polymorphic_alleles_vec = simulation_struct{compute_mode_ctr}.num_simulated_polymorphic_alleles_vec;
                % %                                 simulation_struct{compute_mode_ctr}.old_num_simulated_polymorphic_alleles_vec = save_old_num_simulated_polymorphic_alleles_vec{compute_mode_ctr};
                % %                             end
                % %                         end % function existed before !!!!
                % %
                %                         [het_struct{s_ctr}.plot_x_vec het_struct{s_ctr}.plot_y_vec het_struct{s_ctr}.legend_vec] = ...
                %                             FisherWrightPlotResults(freq_struct, absorption_struct, simulation_struct, ...
                %                             N_vec, expansion_factor, s, mu, num_bins, init_str{1}, iters, ...
                %                             fisher_wright_output_dir, fisher_wright_output_file, tmp_dir, all_s_output_file, mathematica_flag); % plot result
                
                %                         if( (s == s_vec(end)) && (~mathematica_flag) ) % plot different s values together
                %                             my_mkdir( dir_from_file_name(all_s_output_file, 1));
                %                             if(exist(all_s_output_file, 'file'))
                %                                 save(all_s_output_file, 'het_struct', '-append');
                %                             else
                %                                 save(all_s_output_file, 'het_struct');
                %                             end
                %                             plot_inds = [3 4 5];
                %                             FisherWrightPlotResults2(het_struct, s_vec, fisher_wright_output_dir, tmp_dir)
                %                         end
                % %                 end
            end % one plot flag
            ctr = ctr+1; % ctr of init str
        end % loop on init_str
    end % loop on expansion model
    if(unite_results_flag) %
        all_s{s_ctr} = load([fisher_wright_output_file '.mat']);
        all_s{s_ctr}.s = s;
        if(s == s_vec(end)) % reach last s (notice which s_vec element we look at)
            save(all_s_output_file, 'all_s');
        end
        % Plot results
    end % if unite results flag
    s_ctr=s_ctr+1;
end % loop on selection coefficients

for www=0:1
    plot_params.weighted = www
    plot_params.sfs_ind=1; plot_allele_freq(s_vec, {D}, plot_params); % plot all results together for simulation
    plot_params.sfs_ind=2; plot_allele_freq(s_vec, {D}, plot_params);  % same for markov-chain
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test simulation/numerics only for CONSTANT population size where we have also analytic solution
if(test_absorption_time || test_moments)
    N = 500; % take moderate value to let all alleles die
    mu = mu_per_site * (10000 / N); % mutation rate (per nucleotide per generation)
    s=0; % -0.0001;
    expansion_factor = 1.0263; % 1.01; % No EXPANSION! Constant population size
    model_name = 'equilibrium'; % 'expansion 1.01'
    init_str = []; init_str{1} = 'equilibrium'; % CHAGE !!! 'newly_born';
    num_generations = 100; % number of generations to carry with simulations
    iters = 6000; % how many alleles to keep at the end
    compute_mode = 'simulation';
    num_bins = 2*N+1;
    two_side_flag = 0; % Derived allele frequency
    
    % New! create Demography structure
    D.init_pop_size = [N -1]; %   N*2];
    D.generations = [num_generations 1*num_generations]; %  num_generations]; % [num_generations 10*num_generations];
    D.expan_rate = [expansion_factor 1]; %  0.99]; % [expansion_factor 1.0];
    D.index = 1; % only one demography
    D.add_new_alleles = 0; % use old simulation (track only alleles at start) - to be changed !!
    
    N_vec = demographic_parameters_to_n_vec(D, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(test_absorption_time) % Test simulation/numerics only for CONSTANT population size where we have also analytic solution
    [freq_struct_simulation, absorption_struct, simulation_struct, N_vec, simulation_time] = ...
        FisherWrightSimulation([], D, mu, s, init_str{1}, iters, 'simulation', num_bins); % run simulation
    [freq_struct_numeric, absorption_struct_numeric, ~, ~, numeric_time] = ...
        FisherWrightSimulation([], D, mu, s, init_str{1}, iters, 'numeric', num_bins); % run numeric computation
    [freq_struct_moments, absorption_struct_moments, ~, ~, moments_time] = ...
        FisherWrightSimulation([], D, mu, s, init_str{1}, iters, 'moment', num_bins); % run analytic computation using moments
    
    if(debug_figures)
        D_equil = D; D_equil.init_pop_size = D.init_pop_size(1); D_equil.generations = 1000; D_equil.expan_rate = 1;
        N_vec_equil = demographic_parameters_to_n_vec(D_equil, 1);
        [freq_struct_simulation_equil, absorption_struct_equil, simulation_struct_equil] = ...
            FisherWrightSimulation([], D_equil, mu, s, init_str{1}, iters, 'simulation', num_bins); % run simulation
        [freq_struct_numeric_equil, absorption_struct_numeric_equil] = ...
            FisherWrightSimulation([], D_equil, mu, s, init_str{1}, iters, 'numeric', num_bins); % run numeric computation
        [freq_struct_moments_equil, absorption_struct_moments_equil] = ...
            FisherWrightSimulation([], D_equil, mu, s, init_str{1}, iters, 'moment', num_bins); % run analytic computation using moments
        D_equil_end = D_equil; D_equil_end.init_pop_size = N_vec(end-1); D_equil_end.generations = 1000;
        [freq_struct_simulation_equil_end, absorption_struct_equil_end, simulation_struct_equil_end] = ...
            FisherWrightSimulation([], D_equil_end, mu, s, init_str{1}, iters, 'simulation', num_bins); % run simulation
        [freq_struct_numeric_equil_end, absorption_struct_numeric_equil_end] = ...
            FisherWrightSimulation([], D_equil_end, mu, s, init_str{1}, iters, 'numeric', num_bins); % run numeric computation
        [freq_struct_moments_equil_end, absorption_struct_moments_equil_end] = ...
            FisherWrightSimulation([], D_equil_end, mu, s, init_str{1}, iters, 'moment', num_bins); % run analytic computation using moments
                
        D_equil_2000 = D_equil_end; D_equil_2000.init_pop_size = 2000; D_equil_2000.generations = 1000;
        [freq_struct_simulation_equil_2000, absorption_struct_equil_2000, simulation_struct_equil_2000] = ...
            FisherWrightSimulation([], D_equil_2000, mu, s, init_str{1}, iters, 'simulation', num_bins); % run simulation
        [freq_struct_numeric_equil_2000, absorption_struct_numeric_equil_2000, ~, ~, numeric_time_equil_end] = ...
            FisherWrightSimulation([], D_equil_2000, mu, s, init_str{1}, iters, 'numeric', num_bins); % run numeric computation
        [freq_struct_moments_equil_2000, absorption_struct_moments_equil_2000] = ...
            FisherWrightSimulation([], D_equil_2000, mu, s, init_str{1}, iters, 'moment', num_bins); % run analytic computation using moments
    end % if debug_figures
    
    old_cumulative = 0;
    if(old_cumulative)    % take from all generations
        all_p_vec_simulation = accumarray( cell2vec(freq_struct_simulation.x_vec)'+1, cell2vec(freq_struct_simulation.p_vec)'); % average contribution from all generations (why? get the time spent at each allele frequency)
        all_p_vec_numeric = accumarray( cell2vec(freq_struct_numeric.x_vec)'+1, cell2vec(freq_struct_numeric.p_vec)');
        all_p_vec_moments = accumarray( cell2vec(freq_struct_moments.x_vec)'+1, cell2vec(freq_struct_moments.p_vec)'); % NEW! compute #generations from moments! need to debug here!
    else % take last generation
        all_p_vec_simulation = zeros(2*N_vec(end)+1,1);
        all_p_vec_simulation(1:max(freq_struct_simulation.x_vec{end-1})+1) = ...
            accumarray(freq_struct_simulation.x_vec{end-1}'+1, ...
            freq_struct_simulation.p_vec{end-1}') .* ...
            freq_struct_simulation.prob_site_polymorphic_at_end(end);
        all_p_vec_numeric = zeros(2*N_vec(end)+1,1);
        all_p_vec_numeric(1:max(freq_struct_numeric.x_vec{end-1})+1) = ...
            freq_struct_numeric.p_vec{end-1} .* ...
            freq_struct_numeric.prob_site_polymorphic_at_end(end); % This includes the monomorphic state!!!
        all_p_vec_moments = zeros(2*N_vec(end)+1,1);
        all_p_vec_moments(1:max(freq_struct_moments.x_vec{end-1})+1) = ...
            accumarray(freq_struct_moments.x_vec{end-1}'+1, ...
            freq_struct_moments.p_vec{end-1}') * ...
            new_p_poly(end);
        freq_struct_moments.prob_site_polymorphic_at_end(end); % NEW! compute #generations from moments! need to debug here!
        all_p_vec_analytic = 2.*mu.*exp( allele_freq_spectrum((0:2*N) ./ (2*N), s, N, two_side_flag, 'log') ); % compute analytic approxiamtion (valid only for constant population size)
        all_p_vec_analytic_end = 2.*mu.*exp( allele_freq_spectrum((0:2*N_vec(end)) ./ (2*N_vec(end)), s, N_vec(end), two_side_flag, 'log') ); % compute analytic approxiamtion (valid only for constant population size)
        
    end
    
    x_vec = (1:(2*N-1)) ./ (2*N); % take only polymorphic frequencies !!
    x_vec_final = (1:(2*N_vec(end)-1)) ./ (2*N_vec(end)); % change to end-1?
    for cum_flag = 0:1
        for plot_flag = 2:2
            if(plot_flag == 2) % normalize
                bin_size = x_vec(2)-x_vec(1);
                bin_size_final = x_vec_final(2)-x_vec_final(1);
                density_str = '(density)';
            else
                bin_size = 1; bin_size_final = 1;
                density_str = '';
            end
            if(~cum_flag) % Plot PROBABILITY at each allele frequency !!!
                figure; loglog( x_vec_final, all_p_vec_simulation(2:end-1) ./ bin_size_final, 'linewidth', 2 ); hold on;
                loglog( x_vec_final, all_p_vec_numeric(2:end-1) ./ bin_size_final, 'm', 'linewidth', 2 ); % Why divide by # simulations here?? problem with Normalization here!!!
                loglog( x_vec, all_p_vec_analytic(2:end-1) ./ bin_size, 'r' , 'linewidth', 2); % HERE WE MULTIPLY BY FACTOR 2 !!!
                loglog( x_vec_final, all_p_vec_analytic_end(2:end-1) ./ bin_size_final, 'r--' , 'linewidth', 2 );
                loglog( x_vec_final, all_p_vec_moments(2:end-1) ./ bin_size_final, 'c' , 'linewidth', 2 );        % New! add moments based calculations
            else        % cumulative
                figure; semilogx( x_vec_final, cumsum(all_p_vec_simulation(2:end-1)), 'linewidth', 2 ); hold on;
                semilogx( x_vec_final, cumsum(all_p_vec_numeric(2:end-1)), 'm', 'linewidth', 2 ); % Why divide by # simulations here?? problem with Normalization here!!!
                semilogx( x_vec, cumsum(all_p_vec_analytic(2:end-1)), 'r' , 'linewidth', 2); % HERE WE MULTIPLY BY FACTOR 2 !!!
                semilogx( x_vec_final, cumsum(all_p_vec_analytic_end(2:end-1)), 'r--' , 'linewidth', 2 );
                semilogx( x_vec_final, cumsum(all_p_vec_moments(2:end-1)), 'c' , 'linewidth', 2 );        % New! add moments based calculations
                
                
            end
            
            legend({'simulation', 'numeric', 'diffusion-approximation (start)', 'diffusion-approximation (end)', 'Moments'}, ...
                'location', 'southwest', 'fontsize', 14); legend('boxoff');
            xlabel('Derived Allele Frequency'); ylabel(['Prob. ' density_str]);
            title_str = ['N=' num2str(N*1) '->' num2str(N_vec(end-1)) ...
                ', s=' num2str(s,3)  ...
                ', iters ' num2str(simulation_struct.num_simulated_polymorphic_alleles_vec(1)) '->' num2str(iters(1)) ', ' model_name];
            title(str2title([' Prob. ' density_str ' at each allele freq. ' title_str]));
            my_saveas(gcf, fullfile(fisher_wright_output_dir, ...
                ['mean_time_at_each_allele_freq_simulation_vs_diffusion_approx_' density_str(2:end-1)]), 'pdf');
        end
    end
    num_moments = 6; % NEW! plot moments
    analytic_moment_mat = zeros(num_moments, 1); analytic_het_moment_mat = analytic_moment_mat;
    for k=1:num_moments
        analytic_moment_mat(k) = absorption_time_by_selection(s, 1, N, 0, 1, k); % -k-1); % compute (k)-th moment of SFS function (START FROM FIRST MOMENT!!)
        analytic_het_moment_mat(k) = absorption_time_by_selection(s, 1, N, 0, 1, -k-1); % compute (k-1)-th moment of heterozygosity function (START FROM ZERO-TH MOMENT !!!)
    end
    TotalPolymorphicTime = absorption_time_by_selection(s, 1, N, 1/(2*N), 1, 0);
    analytic_moment_mat = analytic_moment_mat ./ TotalPolymorphicTime;
    analytic_het_moment_mat = analytic_het_moment_mat ./ TotalPolymorphicTime;
    combined_moment_mat = [freq_struct_simulation.moments_mat(end-1,:); ...
        freq_struct_numeric.moments_mat(end-1,:); ...
        freq_struct_moments.moments_mat(end-1,:); ...
        analytic_moment_mat'];
    figure; bar(combined_moment_mat'); xlabel('Moment order'); ylabel('Moment value'); title(['Moments of SFS ' title_str]);
    legend({'simulation', 'numeric', 'Moments', 'diffusion-approximation'});   %  plot(freq_struct_simulation.moments_mat(end-1,:)
    combined_het_moment_mat = [freq_struct_simulation.het_moments_mat(end-1,:); ...
        freq_struct_numeric.het_moments_mat(end-1,:); ...
        freq_struct_moments.het_moments_mat(end-1,:); ...
        2*analytic_het_moment_mat'];  % add factor 2 here for: 2x(1-x) vs. just x(1-x)
    figure; bar(combined_het_moment_mat'); xlabel('Het. Moment order'); ylabel('Het. Moment value'); title(['Moments of Het. SFS ' title_str]);
    legend({'simulation', 'numeric', 'Moments', 'diffusion-approximation'});   %  plot(freq_struct_simulation.moments_mat(end-1,:)
    
    if(debug_figures)
        figure; hold on;
        prob_site_polymorphic_at_equilibrium_end = (2*N_vec(end-1)*mu) * 2 * ...
            absorption_time_by_selection(s, 1, N_vec(end-1), 1/(2*N_vec(end-1)), 0.999999999, 0);  % NEW! Factor of two here! fraction of poylmporphic sites at start
        plot(freq_struct_simulation.prob_site_polymorphic_at_equilibrium * ...
            simulation_struct.num_simulated_polymorphic_alleles_vec ./ ...
            simulation_struct.num_simulated_polymorphic_alleles_vec(1));
        plot(absorption_struct_numeric.frac_polymorphic_vec, 'g');
        plot(repmat(prob_site_polymorphic_at_equilibrium_end, length(N_vec)-1, 1), '*r');
        plot(freq_struct_simulation_equil_end.prob_site_polymorphic_at_equilibrium * ...
            simulation_struct_equil_end.num_simulated_polymorphic_alleles_vec ./ ...
            simulation_struct_equil_end.num_simulated_polymorphic_alleles_vec(1), 'b*');
        plot(freq_struct_simulation_equil.prob_site_polymorphic_at_equilibrium * ...
            simulation_struct_equil.num_simulated_polymorphic_alleles_vec ./ ...
            simulation_struct_equil.num_simulated_polymorphic_alleles_vec(1), 'b--');
        plot(absorption_struct_numeric_equil.frac_polymorphic_vec, 'g--'); % run at constant size: old equilibrium
        plot(absorption_struct_numeric_equil_end.frac_polymorphic_vec, 'r--');  % run at constant size: new equilibrium
        
        figure; hold on;
        plot(freq_struct_simulation_equil.prob_site_polymorphic_at_equilibrium * ...
            simulation_struct_equil.num_simulated_polymorphic_alleles_vec ./ ...
            simulation_struct_equil.num_simulated_polymorphic_alleles_vec(1), 'b', 'linewidth', 2);
        plot(absorption_struct_numeric_equil.frac_polymorphic_vec, 'g', 'linewidth', 2); % run at constant size: old equilibrium
        title('N=500')
        figure; hold on;
        plot(freq_struct_simulation_equil_2000.prob_site_polymorphic_at_equilibrium * ...
            simulation_struct_equil_2000.num_simulated_polymorphic_alleles_vec ./ ...
            simulation_struct_equil_2000.num_simulated_polymorphic_alleles_vec(1), 'b', 'linewidth', 2);
        plot(absorption_struct_numeric_equil_2000.frac_polymorphic_vec, 'g', 'linewidth', 2); % run at constant size: old equilibrium
        title('N=2000')
        xlabel('Generation'); ylabel('P_{poly}'); legend('simulation', 'numeric', 'new-equilibrium');
        
        
        figure; hold on; plot(freq_struct_simulation_equil.prob_site_polymorphic_at_end);
        plot(freq_struct_numeric_equil.prob_site_polymorphic_at_end, 'r');
        plot(freq_struct_moments_equil.prob_site_polymorphic_at_end, 'g');
        xlabel('Generation'); ylabel('P_{poly}'); legend('simulation', 'numeric', 'new-equilibrium'); title('equilibrium');
        
        % Plot Prob. (polymorphic) for expansion
        figure; hold on; plot(freq_struct_simulation.prob_site_polymorphic_at_end);
        plot(freq_struct_numeric.prob_site_polymorphic_at_end, 'r');
        plot(freq_struct_moments.prob_site_polymorphic_at_end, 'g');
        xlabel('Generation'); ylabel('P_{poly}'); legend('simulation', 'numeric', 'moments'); title('expansion');
        new_p_poly = 4*2*N*mu*freq_struct_moments.mu_vec_analytic(1,1:end-1) ./ freq_struct_moments.het_moments_mat(:,1)'
        plot(new_p_poly, 'g*');
        
        
        figure; hold on; % plot zero'th moment for equilibrium
        plot(freq_struct_simulation_equil.het_moments_mat(:,1), 'b');
        plot(freq_struct_numeric_equil.het_moments_mat(:,1), 'r');
        plot(freq_struct_moments_equil.het_moments_mat(:,1), 'g');
        plot(freq_struct_moments_equil.mu_vec_analytic(1,1:end-1)' ./ (freq_struct_moments_equil.prob_site_polymorphic_at_end.*N_vec_equil(1:end-1)), 'g--');
        xlabel('Generation'); ylabel('\mu_0'); legend('simulation', 'numeric', 'moments-fitted', 'moments');title('equilibrium');
        
        
        
        figure; hold on; % plot zero'th moment for expansion
        for j=1:num_moments
            subplot(2,3, j); hold on;
            plot(0.5*freq_struct_simulation.het_moments_mat(:,j).*freq_struct_simulation.prob_site_polymorphic_at_end, 'b') ;
            plot(0.5*freq_struct_numeric.het_moments_mat(:,j).*freq_struct_numeric.prob_site_polymorphic_at_end, 'r');
            plot(0.5*freq_struct_moments.het_moments_mat(:,j).*freq_struct_moments.prob_site_polymorphic_at_end, 'g+');
            plot(freq_struct_moments.mu_vec_analytic(j,1:end-1)' .* 4*N*mu, 'co');  % Exacy zero'th moment
            xlabel('Generation'); ylabel('\mu_0'); legend('simulation', 'numeric', 'moments-fitted', 'moments'); title('expansion');
        end
        
        figure; hold on;
        plot(freq_struct_simulation.het_moments_mat(:,2)./freq_struct_simulation.het_moments_mat(:,1), 'b') ;
        plot(freq_struct_numeric.het_moments_mat(:,2)./freq_struct_numeric.het_moments_mat(:,1), 'r') ;
        plot(freq_struct_moments.het_moments_mat(:,2)./freq_struct_moments.het_moments_mat(:,1), 'g') ;
        plot(freq_struct_moments.mu_vec_analytic(2,:) ./ freq_struct_moments.mu_vec_analytic(1,:), 'g*');
        xlabel('Generation'); ylabel('\mu_0'); legend('simulation', 'numeric', 'moments-fitted'); title('expansion'); title('\mu_1 / \mu_0 ratio');
        
        figure; hold on;
        plot(freq_struct_simulation.het_moments_mat(:,1), 'b') ;
        plot(freq_struct_numeric.het_moments_mat(:,1), 'r') ;
        plot(freq_struct_moments.het_moments_mat(:,1), 'g') ;
        plot(freq_struct_moments.mu_vec_analytic(1,1:end-1)  .* 8*N*mu ./ freq_struct_moments.prob_site_polymorphic_at_end', 'g*') ;
        xlabel('Generation'); ylabel('\mu_0'); legend('simulation', 'numeric', 'moments-fitted'); title('expansion'); title('\mu_0 ');
        
        
        
        figure;
        semilogx(freq_struct_numeric.x_vec{1}, cumsum(freq_struct_numeric.p_vec{1})); hold on;
        semilogx(freq_struct_simulation.x_vec{1}, cumsum(freq_struct_simulation.p_vec{1}), 'g');
        semilogx(freq_struct_moments.x_vec{1}, cumsum(freq_struct_moments.p_vec{1}), 'r');
        xlabel('x'); ylabel('\Phi(x)'); legend('numeric',  'simulation', 'moments'); title('cumulative \Phi at t=1');
        
        figure;
        semilogx(freq_struct_numeric.x_vec{end-1}, cumsum(freq_struct_numeric.p_vec{end-1}));  hold on;
        semilogx(freq_struct_simulation.x_vec{end-1}, cumsum(freq_struct_simulation.p_vec{end-1}), 'g');
        semilogx(freq_struct_moments.x_vec{end-1}, cumsum(freq_struct_moments.p_vec{end-1}), 'r');
        xlabel('x'); ylabel('\Phi(x)'); legend('numeric', 'simulation', 'moments'); title(['cumulative \Phi at t=' length(freq_struct_numeric.x_vec)-1]);
        
        freq_struct_numeric.het_moments_mat(1,:)
        freq_struct_moments.het_moments_mat(1,:)
        
        freq_struct_numeric.het_moments_mat(end-1,:)
        freq_struct_moments.het_moments_mat(end-1,:)
        
        num = load('numeric_file.mat'); mom = load('moments_file.mat');
        figure; hold on;
        plot(num.new_vec); plot(num.old_vec, 'r');
        plot(mom.new_vec, 'b--'); plot(mom.old_vec, 'r--');
        legend('new-num', 'old-num', 'new-mom', 'old-mom');
        
        
        
    end
    
    % % %     % Move to continuous limit - plot densities:
    % % %     legend('simulation', 'diffusion-approximation', 'half-diffusion-approximation');
    % % %     xlabel('Derived Allele Frequency'); ylabel('Num. Generations Spent (Density)');
    % % %     title_str = ['N=' num2str(N*1) ...
    % % %         ', s=' num2str(s,3)  ...
    % % %         ', iters ' num2str(simulation_struct.num_simulated_polymorphic_alleles_vec(1)) '->' num2str(iters(1))];
    % % %     title(str2title([' Mean # Generations (density) spent at each allele freq. ' title_str]));
    % % %     my_saveas(gcf, 'mean_time_at_each_allele_freq_simulation_vs_diffusion_approx', 'pdf');
    
    %    bar_mat = [all_p_vec_simulation(2:end-1) freq_struct_numeric.p_vec{end-1}(2:(end-1))']
    %    figure; hold on; bar(bar_mat); legend('simulation', 'numeric');
    %    figure; hold on; plot(all_p_vec_simulation(2:end-1), freq_struct_numeric.p_vec{end-1}(2:(end-1)), '*'); xlabel('simulation'); ylabel('numeric');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute moments of SFS
if(test_moments)
    max_k = 5; % how many moments to compute
    D.expan_rate = 1.1;
    D.index = 1;
    D.generations = 100;
    N_vec = demographic_parameters_to_n_vec(D, 1);
    [mu_vec_expansion_analytic, mu_vec_equilibrium] = FisherWright_Compute_SFS_Moments(N_vec, 0, max_k); % compute moments with Formulas from Ewens
    
    % Estimate density using the max-entropy method:
    x_vec = (1:(2*N-1)) ./ (2*N);
    [lambda_max_ent, g_het_max_ent, entr_max_ent] = ...  % Fit density using moments % , lambda0)
        me_dens2(mu_vec_expansion_analytic(2:end) ./ mu_vec_expansion_analytic(1), x_vec); % normalize by 0th moment
    
    f_max_ent = zeros(1, 2*N-1);
    for k=1:(max_k)
        f_max_ent = f_max_ent + lambda_max_ent(k) .* x_vec .^ (k-1);
    end
    f_max_ent = exp(-f_max_ent) ./ (x_vec .* (1-x_vec)); % Compute f. How to normalize?
    figure; plot(x_vec, f_max_ent); xlabel('f'); ylabel('phi_0(f) density');
    
    iters = 2000; % Compute moments using simulation
    [freq_struct_simulation, absorption_struct, simulation_struct, N_vec, simulation_time] = ...
        FisherWrightSimulation([], D, mu, s, init_str{1}, iters, 'simulation', num_bins); % run simulation
    final_x_vec_simulations = freq_struct_simulation.x_vec{num_generations-1} ./ (2*N)
    % Compute moments
    mu_vec_expansion_simulations = zeros(max_k, 1);
    for k=1:max_k
        mu_vec_expansion_simulations(k) = moment_hist(final_x_vec_simulations, ...
            final_x_vec_simulations .* (1-final_x_vec_simulations) .* freq_struct_simulation.p_vec{num_generations-1}, k, 0);
    end
    
    
    % Compare two:
    
end

total_time = cputime - total_time

% % Check binomial and Gaussian approximations
% n = 10000; p = 0.005; M = 10000;
% y_binom = binornd(n, p, [M, 1]);
% y_poiss = poissrnd(n*p, [M, 1]);
% y_normal = normrnd(n*p, sqrt(n*p*(1-p)), [M, 1]);
%
% figure; hold on;
% plot(sort(y_binom), [1:M]./M);
% plot(sort(y_poiss), [1:M]./M, 'g');
% plot(sort(y_normal), [1:M]./M, 'r');
% legend('binom', 'poisson', 'normal');
%








