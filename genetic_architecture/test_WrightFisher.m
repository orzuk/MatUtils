% Script for testing calculations with the Wright-Fisher model
run_sim_flag = 0; % compute site frequency distribution for each s
run_two_expansions_flag = 0; % a mode within simulation, running a two-stage expansion model
one_plot_flag = 0; % make one plot for each s
unite_results_flag = 0; % make one plot for all s values together
save_in_mathematica = 0; % save in mathematica format for Eric
test_absorption_time = 1; % compare absorption time for simulations vs. analytical formula


AssignGeneralConstants;
N = 500; % population size
num_generations_vec = [200 10]; % two-stage model: slow and fast expansion % 2500; % 50;
expansion_factor_vec = [1.005 1.05]; % 1.1; % 1.02; % two-stage moe: growth in population size


% s_vec = -[0 logspace(-6, -1, 11)]; % take log-space
s_vec = -[0 0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01]; %  0.05 0.1]; % -0.00000001; % selection coefficient
compute_mode = 'simulation'; % 'simulation'; % 'simulation';  % 'simulation'; % 'numeric'; % 'simulation'; % 'numeric'; % how to advance calculation
mu = 2*10^(-8) * (10000 / N); % mutation rate (per nucleotide per generation)
init_str = 'equilibrium'; % 'equilibrium' 'newly_born'; % start at newly born allele or equilibrium
iters = [20000 20]; % relevant only for simulation. Spend more iterations on equilibrium (new alleles take more time per iteration)

if(~save_in_mathematica)
    num_bins = 2000; % for plotting histograms
    mathematica_str = '';
    mathematica_flag = 0;
else
    num_bins = [0 logspace(-6, 0, 601)]; % use logarithmic bins  (for Eric)
    mathematica_str = '_mathematica';
    mathematica_flag = 1;
end
fisher_wright_output_dir = '../../common_disease_model/docs/pnas/power_paper/figs/SFS';


s_ctr = 1; % counter of selection coefficient
het_struct = cell(length(s_vec), 1);

total_time = cputime;
for s = s_vec(1) %   1:end-2) % (2:end) % s_vec(end) % s_vec(3:end) % (end-1) % s_vec % loop on different selection coefficients
    close all; run_s = s
    frac_polymorphic = 2*N*mu * absorption_time_by_selection(s, 1, N, 1/(2*N), 0.999999999, 0);
    for expansion_model = 0:0 % run_two_expansions_flag % allow one rate and two-rate exponential expansions
        ctr=1; % counter on init str
        
        if(expansion_model) % here run two-stage expansion
            expansion_str = '_two-stage-expansion';
            num_generations = sum(num_generations_vec); % simulate a few stages together
            expansion_factor = zeros(num_generations, 1);
            expansion_ctr = 1;
            for i=1:length(expansion_factor_vec)
                expansion_factor(expansion_ctr:(expansion_ctr+num_generations_vec(i)-1)) = expansion_factor_vec(i);
                expansion_ctr = expansion_ctr + num_generations_vec(i);
            end % loop on different expansions
            tmp_dir = (['N_' num2str(N) '_two_stage_expansion_k_' num2str(num_generations) '_s_' num2str(s,3)]);
            
        else % here run one stage exponential expansion
            expansion_str = '';
            expansion_factor = expansion_factor_vec(1);
            num_generations = num_generations_vec(1);
            tmp_dir = (['N_' num2str(N) '_expansion_' num2str(expansion_factor,3) '_k_' num2str(num_generations) '_s_' num2str(s,3)]);
        end
        tmp_dir = strrep(tmp_dir, '.', '_');
        for init_str = {'equilibrium', 'newly_born'} % , 'equilibrium' 'newly_born'} % list two distributions separately
            fisher_wright_output_file = fullfile(fisher_wright_output_dir, tmp_dir, ['fisher_wright_expansion_' init_str{1}]);
            all_s_output_file = fullfile(fisher_wright_output_dir, ['fisher_wright_expansion_all_s_' ...
                init_str{1} expansion_str mathematica_str]);
            
            if(run_sim_flag)
                compute_mode_ctr=1; freq_struct = cell(2,1); absorption_struct = cell(2,1); simulation_struct = cell(2,1); simulation_time = zeros(2,1)
                for compute_mode = {'simulation', 'numeric'}
                    [freq_struct{compute_mode_ctr} absorption_struct{compute_mode_ctr} simulation_struct{compute_mode_ctr} ...
                        N_vec{compute_mode_ctr} simulation_time(compute_mode_ctr)] = ...
                        FisherWrightSimulation(N, mu, s, num_generations, expansion_factor, init_str{1}, iters(ctr), compute_mode{1}, num_bins); % run simulation
                    compute_mode_ctr=compute_mode_ctr+1;
                end
                
                my_mkdir(fullfile(fisher_wright_output_dir, tmp_dir));
                save([fisher_wright_output_file '.mat'], ...
                    'freq_struct', 'absorption_struct', 'simulation_struct', 'N_vec', 'simulation_time');
            else
                if(one_plot_flag) % no need to load if not plotting anything
                    load([fisher_wright_output_file '.mat']);
                end
            end
            if(one_plot_flag) % one plot
                switch init_str{1}
                    case 'equilibrium' % This assumes equilibrium is called first !!!
                        for compute_mode_ctr=1:2
                            save_all_old_x_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.x_vec{num_generations}; %  .* (2*N_vec(num_generations));
                            save_all_old_p_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.p_vec{num_generations};
                            save_all_old_het_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.het_vec{num_generations};
                            
                            save_old_total_het_at_each_generation_vec{compute_mode_ctr} = freq_struct{compute_mode_ctr}.total_het_at_each_generation_vec
                            if(compute_mode_ctr==1) % simulations
                                save_old_num_simulated_polymorphic_alleles_vec{compute_mode_ctr} = simulation_struct{compute_mode_ctr}.num_simulated_polymorphic_alleles_vec;
                            end
                        end
                    case 'newly_born'
                        for compute_mode_ctr=1:2
                            freq_struct{compute_mode_ctr}.all_old_x_vec = save_all_old_x_vec{compute_mode_ctr};
                            freq_struct{compute_mode_ctr}.all_old_p_vec = save_all_old_p_vec{compute_mode_ctr};
                            freq_struct{compute_mode_ctr}.all_old_het_vec = save_all_old_het_vec{compute_mode_ctr};
                            freq_struct{compute_mode_ctr}.old_total_het_at_each_generation_vec = save_old_total_het_at_each_generation_vec{compute_mode_ctr};
                            
                            if(compute_mode_ctr==1) % simulations
                                simulation_struct{compute_mode_ctr}.new_num_simulated_polymorphic_alleles_vec = simulation_struct{compute_mode_ctr}.num_simulated_polymorphic_alleles_vec;
                                simulation_struct{compute_mode_ctr}.old_num_simulated_polymorphic_alleles_vec = save_old_num_simulated_polymorphic_alleles_vec{compute_mode_ctr};
                            end
                        end
                        [het_struct{s_ctr}.plot_x_vec het_struct{s_ctr}.plot_y_vec het_struct{s_ctr}.legend_vec] = ...
                            FisherWrightPlotResults(freq_struct, absorption_struct, simulation_struct, ...
                            N_vec, expansion_factor, s, mu, num_bins, init_str{1}, iters, ...
                            fisher_wright_output_dir, fisher_wright_output_file, tmp_dir, all_s_output_file, mathematica_flag); % plot result
                        
                        
                        if( (s == s_vec(end)) && (~mathematica_flag) ) % plot different s values together
                            my_mkdir( dir_from_file_name(all_s_output_file, 1));
                            if(exist(all_s_output_file, 'file'))
                                save( [all_s_output_file], 'het_struct', '-append');
                            else
                                save( [all_s_output_file], 'het_struct');
                            end
                            plot_inds = [3 4 5];
                            
                            
                            FisherWrightPlotResults2(het_struct, s_vec, fisher_wright_output_dir, tmp_dir)
                            
                            
                        end
                end
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



if(test_absorption_time) % Test simulation/numerics only for CONSTANT population size where we have also analytic solution
    N = 500; % take moderate value to let all alleles die 
    mu = 2*10^(-8) * (10000 / N); % mutation rate (per nucleotide per generation)
    s=0; % -0.0001;
    expansion_factor = 1; % No EXPANSION! Constant population size
    init_str{1} = 'newly_born';
    num_generations = 800; % number of generations to carry with simulations 
    iters = 30; % how many alleles to keep at the end
    compute_mode = 'simulation';
    num_bins = 2*N+1;
    two_side_flag = 0; % Derived allele frequency
    %     [q x_vec p_vec p_vec_analytic het_vec N_vec ...
    %         absorption_time_vec fixation_time_vec loss_time_vec total_het_vec ...
    %         frac_polymorphic_vec prob_fixation frac_alleles_kept_vec frac_het_kept_vec prob_site_polymorphic ...
    %         all_new_x_vec all_new_p_vec all_new_het_vec num_simulated_polymorphic_alleles_vec simulation_time] = ...
    [freq_struct_simulation absorption_struct simulation_struct N_vec simulation_time] = ...
        FisherWrightSimulation(N, mu, s, num_generations, expansion_factor, init_str{1}, iters, 'simulation', num_bins); % run simulation
    [freq_struct_numeric absorption_struct_numeric simulation_struct_numeric N_vec numeric_time] = ...
        FisherWrightSimulation(N, mu, s, num_generations, expansion_factor, init_str{1}, iters, 'numeric', num_bins); % run numeric computation
        
    all_p_vec_simulation = accumarray( cell2vec(freq_struct_simulation.x_vec)'+1, cell2vec(freq_struct_simulation.p_vec)'); % average contribution from all generations (why? get the time spent at each allele frequency)
    all_p_vec_numeric = accumarray( cell2vec(freq_struct_numeric.x_vec)'+1, cell2vec(freq_struct_numeric.p_vec)');
    x_vec = (1:(2*N-1)) ./ (2*N);
        
    for plot_flag = 1:2
        if(plot_flag == 2) % normalize 
            bin_size = x_vec(2)-x_vec(1);
            density_str = '(density)';
        else
            bin_size = 1;
            density_str = '';
        end
        figure; loglog( x_vec, (all_p_vec_simulation(2:end-1) ./ 1) ./ 1 ); hold on;
%        figure; loglog( x_vec, (all_p_vec_simulation(2:end-1) ./ simulation_struct.num_simulated_polymorphic_alleles_vec(1)) ./ bin_size ); hold on;
        loglog(x_vec, (all_p_vec_numeric(2:end-1) ./ simulation_struct.num_simulated_polymorphic_alleles_vec(1)) ./ bin_size, 'm' ); % Why divide by # simulations here?? 
        all_p_vec_analytic = exp( allele_freq_spectrum((0:2*N) ./ (2*N), s, N, two_side_flag, 'log') ); % compute analytic approxiamtion (valid only for constant population size)
        loglog( x_vec, (4*N*mu) * (2.*all_p_vec_analytic(2:end-1) ./ (2*N)) ./ bin_size, 'r' , 'linewidth', 2); % HERE WE ADD FACTOR 2 !!! 
        loglog( x_vec, (4*N*mu) * (all_p_vec_analytic(2:end-1) ./ (2*N)) ./ bin_size, 'g--' , 'linewidth', 2 );
        
        %         figure; loglog( x_vec, all_p_vec(2:end-1) ./ simulation_struct.num_simulated_polymorphic_alleles_vec(1) ); hold on;
        %         all_p_vec_analytic = exp( allele_freq_spectrum((0:2*N) ./ (2*N), s, N, two_side_flag, 'log') ); % compute analytic approxiamtion (valid only for constant population size)
        %         loglog( x_vec, 2.*all_p_vec_analytic(2:end-1) ./ (2*N), 'r' , 'linewidth', 2);
        %         loglog( x_vec, all_p_vec_analytic(2:end-1) ./ (2*N), 'g--' , 'linewidth', 2 );
        legend('simulation', 'numeric', 'diffusion-approximation', 'half-diffusion-approximation');
        xlabel('Derived Allele Frequency'); ylabel(['Num. Generations Spent ' density_str]);
        title_str = ['N=' num2str(N*1) ...
            ', s=' num2str(s,3)  ...
            ', iters ' num2str(simulation_struct.num_simulated_polymorphic_alleles_vec(1)) '->' num2str(iters(1))];
        title(str2title([' Mean # Generations ' density_str ' newly born allele spent at each allele freq. ' title_str]));
        my_saveas(gcf, ['mean_time_at_each_allele_freq_simulation_vs_diffusion_approx_' density_str(2:end-1)], 'pdf');
    end
    
    % % %     % Move to continuous limit - plot densities:
    % % %     legend('simulation', 'diffusion-approximation', 'half-diffusion-approximation');
    % % %     xlabel('Derived Allele Frequency'); ylabel('Num. Generations Spent (Density)');
    % % %     title_str = ['N=' num2str(N*1) ...
    % % %         ', s=' num2str(s,3)  ...
    % % %         ', iters ' num2str(simulation_struct.num_simulated_polymorphic_alleles_vec(1)) '->' num2str(iters(1))];
    % % %     title(str2title([' Mean # Generations (density) spent at each allele freq. ' title_str]));
    % % %     my_saveas(gcf, 'mean_time_at_each_allele_freq_simulation_vs_diffusion_approx', 'pdf');
end


total_time = cputime - total_time
