% Compute allele frequency distribution using Fisher-Wright model with changing population size.
% Input:
% D - structure with demographic models
% s - selection coefficient. Should be NEGATIVE for deleterious alleles
% compute_flag - 'simulation' (default) or 'moments' (computation based on moments)
% n_sample - # of individuals in a SAMPLE (default: = population size at END)
% mu - regional mutation rate (default: mu for one site)
%
% Output:
% x_vec - vector of x values (allele frequencies) at each generation
% p_vec - vector of their frequencies at each generation
% L_correction_factor - correction factor for total mutation rate
% compute_time - total time it took to run
% k_vec - alleles in sample
% n_vec - sample sizes
% weights_vec - weight of each allele (each allele can represent multiple alleles)
%
function [x_vec, p_vec, L_correction_factor, compute_time, k_vec, n_vec, weights_vec] = ...
    compute_allele_freq_spectrum_from_demographic_model(D, s, compute_flag, n_sample, mu)

if(~isfield(D, 'save_flag'))
    D.save_flag = 0; 
end
if(~ischar(compute_flag)) % allow for structure of compute parameters
    smooth_params.smooth = compute_flag.smooth;
    compute_flag = compute_flag.method;
else
    smooth_params = [];
end
smooth_params.knots = 10; smooth_params.plot = 0;

if(~isscalar(s)) % NEW! allow to fit multiple s values using a surface fitting module
    s_vec = s; num_s = length(s_vec);  [x_vec_cell, p_vec_cell, s_vec_cell] = deal(cell(num_s, 1)); 
    for i_s = 1:num_s
        sprintf('Run selection s=%f', s_vec(i_s))
        s=s_vec(i_s)
        [x_vec_cell{i_s}, p_vec_cell{i_s}, L_correction_factor, compute_time] = ...
            compute_allele_freq_spectrum_from_demographic_model(D, s, compute_flag); % , n_sample, mu);
        s_vec_cell{i_s} = repmat(s, 1, length(x_vec_cell{i_s})); 
    end
    save(['temp_surface.' D.name '.mat'], 'x_vec_cell', 's_vec_cell', 'p_vec_cell', 'smooth_params', 'D'); 
    
    temp_DEBUG=1;
    if(temp_DEBUG)
        s_vec = mean_cell(s_vec_cell);
        D2=D; D2.s_grid = s_vec; 
        D2.SFS.x_vec = x_vec_cell;
        D2.SFS.p_vec = p_vec_cell;
        plot_params.figure_type = 1; plot_params.figs_dir = []; plot_params.hist = 1; plot_params.xlim = [10^(-4) 1];
        plot_params.cum=1; plot_params.weighted = 1; plot_params.normalize=1; plot_params.font_size=8; % plot cumulative weighted allele frequency distribution
        plot_allele_freq(s_vec, {D2}, plot_params);
    end     
        
    % Perform smoothing with monotonicity constraints:
    x_vec = unique( [x_vec_cell{:}] ); num_x = length(x_vec); 
    p_mat = zeros(num_s, num_x);
    for i_s = 1:length(s_vec)
        [~, I, J] = intersect(x_vec, x_vec_cell{i_s});
        p_mat(i_s,I) = p_vec_cell{i_s}(J);
    end
    %    p_mat = gridfit(double([x_vec_cell{:}]) ,[s_vec_cell{:}], [p_vec_cell{:}], double(x_vec), -s_vec); % fit surface
    if(isfield(D, 's_grid'))
        smooth_params.y_fit = abs(D.s_grid);
    end
    save(['temp_surface.' D.name '.mat'], '-append', 'x_vec', 's_vec', 'p_mat', 'smooth_params', 'D'); 
%    load('temp_surface.mat'); 
    [x_vec, ~, p_vec] = fit_monotonic_surface(x_vec, abs(s_vec), p_mat, smooth_params);  % constraints
    %	
    compute_time=cputime
    return;
else    
    smooth_params.y_fit = s; % fit one s
end% if s is not scalar 


compute_time=cputime
if(~exist('compute_flag', 'var') || isempty(compute_flag))
    compute_flag = 'simulation';
end
% if(~exist('init_str', 'var') || isempty(init_str)) % default is start at equilibrium
init_str = 'equilibrium';
%end
if(~exist('mu', 'var') || isempty(mu))
    AssignRVASConstants;
    mu = mu_per_site; % set default mutation rate (per-nucleotide per-generation)
end
if(~isfield(D, 'iters'))
    D.iters = 10000; % number of alleles to simulate (start low to save time. As we refine demography fitting we increase this number)
end
D.num_bins = 100; % used for binning in Fisher Right simulation. 100 is too little??? 
D.compute_absorb = 0; % no need for extra computation!!!
N_vec = demographic_parameters_to_n_vec(D, D.index); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation
if(~exist('n_sample', 'var') || isempty(n_sample))
    n_sample =  2*N_vec(end-1); % Get last population size
end
weights_vec = 1;
num_final_generations = length(N_vec)-1; % simulation at the end
L_correction_factor=[];

switch compute_flag
    case {'simulation', 'simulations', 'numeric'}
        max_num_alleles = 20000; % set maximum to save time
        [freq_struct, ~, simulation_struct, N_vec, simulation_time] = ... % New: separate output to different structures
            FisherWrightSimulation([], D, mu, s, init_str, D.iters, compute_flag, D.num_bins);
        [freq_struct_num, ~, simulation_struct_num, N_vec_num, simulation_time_num] = ... % New: separate output to different structures
            FisherWrightSimulation([], D, mu, s, init_str, D.iters, 'numeric', D.num_bins);
        p_stationary = allele_freq_spectrum((1:(2*N_vec(1)-1))./(2*N_vec(1)), s, N_vec(1), 0, 'linear');
        p_stationary_num = allele_freq_spectrum_numeric((1:(2*N_vec(1)-1))./(2*N_vec(1)), s, N_vec(1), 0, 'linear');
        p_stationary_num = p_stationary_num ./ sum(p_stationary_num); % normalize
        fprintf('Fisher-Wright simulation time was %f\n', simulation_time);
        if(D.save_flag) % save input and output to file for debugging 
            save(['debug_SFS_' D.name num2str(D.cond_on_polymorphic_flag) '.mat'], ...
                'D', 'mu', 's', 'init_str', 'compute_flag', 'freq_struct', 'simulation_struct', 'N_vec', 'simulation_time');
        end
        M = FisherWright_ComputeMarkovMatrix(N_vec(1), s, 'exact', 1); % compute Markov matrix 
        M(1,2)=1; M(1,1)=1-1; M(end,2)=1; M(end,end)=1-1; % add small mutation to make process ergodic 

        g=1; % pick generation and plot
        p_vec_sim = zeros(size(p_stationary_num)); p_vec_sim(freq_struct.x_vec{g}(2:(end-1))) = freq_struct.p_vec{g}(2:(end-1))' ./ sum(freq_struct.p_vec{g}(2:(end-1))');
        figure; bar([p_vec_sim' freq_struct_num.p_vec{g}(2:(end-1))' p_stationary_num']); legend('sim', 'numeric', 'stationary'); title(['Generation ' num2str(g) ', s=' num2str(s)]);
        for g=1:length(freq_struct.p_vec)
            mean_sim_vec(g) = mean_hist(freq_struct.x_vec{g}, freq_struct.p_vec{g}); 
        end
        figure; plot(mean_sim_vec); hold on; plot(mean_hist(freq_struct.x_vec{g}(2:end-1), p_stationary_num), '*r'); 
        
%         if(isfield(simulation_struct, 'q'))
%             [sorted_q, sort_perm] = sort(simulation_struct.q(:,end));
%             figure; semilogx(sorted_q, cumsum(simulation_struct.weights(sort_perm)) ./ sum(simulation_struct.weights(sort_perm)));
%         end
        x_vec = freq_struct.x_vec{end-1}; % why don't take last one?
        p_vec = freq_struct.p_vec{end-1}; 
        
        L_correction_factor = simulation_struct.L_correction_factor;
%        [x_vec, ~, p_vec] = fit_monotonic_surface(x_vec, abs(s), p_vec, smooth_params);  % constraints
        p_vec = freq_struct.p_vec{end-1} ./ sum(freq_struct.p_vec{end-1});
        
        if(nargout > 4) % output k_vec, n_vec ...
            if(~isfield(simulation_struct, 'weights'))
                simulation_struct.weights = [];
            end
            % currently: round p_vec to integers! (why? could be innaccurate, and also have many alleles for large N)
            p_vec_counts = round(p_vec .* simulation_struct.num_simulated_polymorphic_alleles_vec(end));
            if(sum(p_vec_counts) > max_num_alleles)
                for M=(sum(p_vec_counts)-max_num_alleles):-1:0 % determine where to chop
                    if(sum(min(max(p_vec_counts)-M, p_vec_counts)) > max_num_alleles)
                        break;
                    end
                end
                weights_vec = p_vec_counts ./ min(p_vec_counts, max(p_vec_counts)-M+1);
                p_vec_counts = min(p_vec_counts, max(p_vec_counts)-M+1);
            else
                weights_vec =1;
            end
            allele_freq_vec = hist_to_vals(x_vec, p_vec_counts); % compute alleles at sample ? why not p_vec? 
            num_alleles = length(allele_freq_vec);
            pop_to_sample_t = cputime;
            k_vec = population_to_sample_allele_freq(allele_freq_vec, 2*N_vec(end-1), n_sample); % simulate a sample from population
            pop_to_sample_t = cputime-pop_to_sample_t; % note: we don't always need k_vec, n_vec !!!
            fprintf('Converted %d alleles to sample freq. time=%f\n', num_alleles, pop_to_sample_t);
            n_vec = repmat(n_sample, num_alleles, 1);
        end % if nargout > 4
    case 'moments'
        N_vec = demographic_parameters_to_n_vec(D, 1);
        mu_vec_expansion_analytic = FisherWright_Compute_SFS_Moments(N_vec, 0, max_k); % compute moments with Formulas from Ewens
        
        % Estimate density using the max-entropy method? NO! just compute moments !
        x_vec = (1:(2*N-1)) ./ (2*N); % use initial pop size for resolution
        lambda_max_ent = ...  % Fit density using moments % , lambda0)
            me_dens2(mu_vec_expansion_analytic(2:end) ./ mu_vec_expansion_analytic(1), x_vec); % normalize by 0th moment
        p_vec = zeros(1, 2*N-1);
        for k=1:(max_k)
            p_vec = p_vec + lambda_max_ent(k) .* x_vec .^ (k-1);
        end
        p_vec = exp(-p_vec) ./ (x_vec .* (1-x_vec)); % Compute f. How to normalize?
end

% % % % New! smooth SFS  (should we do it at sample or at population level?)
% % % if(~isempty(smooth_params))
% % %     [x_vec2, p_vec2] = fit_monotonic_curve(log(x_vec(x_vec>0)), p_vec(x_vec>0), smooth_params);
% % %     x_vec2 = exp(x_vec2);
% % %     
% % %     figure; semilogx(x_vec, p_vec, '*'); hold on; % Temp: plot
% % %     semilogx(x_vec2, p_vec2, 'r');
% % %     legend({'data', 'fitted'});
% % % end

compute_time=cputime-compute_time;



