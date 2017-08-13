% Compute allele frequency distribution using Fisher-Wright model with changing population size.
% Input:
% D - structure with demographic models
% s - selection coefficient. NOTE: Should be NEGATIVE for deleterious alleles
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

if(~isstring(compute_flag)) % allow for structure of compute parameters
	smooth_params = compute_flag.smooth;
	compute_flag = compute_flag.method; 
end

if(~isscalar(s)) % NEW! allow to fit multipole s values using a surface fitting module
	s_vec = s; 
	for i_s = 1:length(s_vec)
		p_mat = zeros(num_x, num_s); 
		s=s_vec(s); 
		[x_vec, p_mat(:,i_s), L_correction_factor, compute_time, k_vec, n_vec, weights_vec] = ...
		    compute_allele_freq_spectrum_from_demographic_model(D, s, compute_flag, n_sample, mu);		
	end
	% Perform smoothing with monotonicity constraints: 
    %	[zgrid,xgrid,ygrid] = gridfit(x_mat,s_mat, p_mat,  x_grid,s_grid); % fit monotonic surface 				
end


compute_time=cputime;
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
if(~isfield('iters', D))
    D.iters = 1000; % number of alleles to simulate (start low to save time. As we refine demography fitting we increase this number)
end
D.num_bins = 100; % used for binning in Fisher Right simulation
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
        fprintf('Fisher-Wright simulation time was %f\n', simulation_time);        
        x_vec = freq_struct.x_vec{end-1}; % why don't take last one?
        p_vec = freq_struct.p_vec{end-1};
        L_correction_factor = simulation_struct.L_correction_factor;
        
        if(nargout > 4) % output k_vec, n_vec ...
            if(~isfield(simulation_struct, 'weights'))
                simulation_struct.weights = [];
            end
            % currently: round p_vec to integers! (could be innaccurate, and also have many alleles for large N)
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
            allele_freq_vec = hist_to_vals(x_vec, p_vec_counts); % compute alleles at sample ?
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

% New! smooth SFS  (should we do it at sample or at population level?)
if(smooth_params)
	slm = slmengine(x_vec, p_vec,'plot','on','knots',10,'decreasing','on'); % 'leftslope',0,'rightslope',0);
	p_vec2 = slmeval(x_vec, slm); 

	% Temp: plot 
end

compute_time=cputime-compute_time;
