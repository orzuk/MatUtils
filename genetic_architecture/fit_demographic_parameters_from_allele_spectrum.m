% Fit a simple demographic model to population allele frequency spectrum
% We assume that observed alleles are not under selection (neutral)
%
% Input:
% k_vec - number of dervied allele carriers observed in population for each allele
% n_vec - number of individuals profiled in population for each allele
% mu - total mutation rate for region (assume this is known for now) 
% 
% Output:
% D - demographic model (expansion, population size etc.)
% max_LL - log-likelihood of data for spectrum
%
function [D, max_LL] = fit_demographic_parameters_from_allele_spectrum(k_vec, n_vec, mu, L_correction_factor, D_opt)

AssignGeneralConstants;    AssignRVASConstants; 

% We fit expansion model from allele-frequency data, assuming that all alleles are neutral
s = 0; % Assume no selection (synonymous)
alpha=0; % No mixture. Everything is neutral (synonymous)
beta=[];  % effect size on phenotype IRRELEVANT! (we use only genpotypes)

X = [vec2row(k_vec) vec2row(n_vec)]'; % Represent counts in a packed form

% Set different demographic models:
D.init_pop_size_vec{1} = 500; % 1. ancestral population
D.init_pop_size_vec{2} = round(logspace(1, 4, 10)); % 2. bottleneck size
D.init_pop_size_vec{3} = -1; % 3. population after bottleneck
%D.num_params_vec(1) = length(D.init_pop_size_vec{2});

D.generations_vec{1} = 500; % 4. burnout time
D.generations_vec{2} = 10;  % 5. number of generations since bottleneck
D.generations_vec{3} = round(logspace(1, 4, 10));  % 6. number of generations for expansion
%D.num_params_vec(2) = length(D.generations_vec{3});

D.expan_rate_vec{1} = 1; % fixed ancestral population size
D.expan_rate_vec{2} = 1; % fixed buttleneck size
D.expan_rate_vec{3} = logspace(0.001, 0.04, 10); % 5. expansion rate per generation
%D.num_params_vec(3) = length(D.expan_rate_vec{3});

D.num_params_vec = [length_cell(D.init_pop_size_vec) length_cell(D.generations_vec) length_cell(D.expan_rate_vec)];
D.num_stages = length(D.init_pop_size_vec);

%D.num_bottleneck_size = length(D.bottleneck_size);
%D.num_bottleneck_size = length(D.bottleneck_size);
%D.num_init_pop_size = length(D.bottleneck_size);
%D.num_expansion_rate = length(D.expansion_rate);


D.num_params = prod(D.num_params_vec); % how many parameters to enumerate
compute_flag = 'simulation'; % numeric 

log_like_mat = zeros(size(D));

filter_by_moments=1; num_moments=1; % First filter by moments !!!
if(filter_by_moments)
    moments_time=cputime;
    het_moment_mat_all_models = zeros(D.num_params, num_moments);
    init_N = zeros(D.num_params, 1);
    for i=1:D.num_params
        if(mod(i, 100) == 0)
            sprintf('Fitting %ld out of %ld demographic model', i, D.num_params)
            cur_moment_time = cputime - moments_time
        end
        D.index = i;
        i_vec = myind2sub(D.num_params_vec, length(D.num_params_vec), i); % create vector of indices
        if(i_vec(2*D.num_stages) == length(D.generations_vec{end})) % get highest number of generations

            j_vec = i_vec;
            for j=length(D.generations_vec{end}):-1:1 % loop on all smaller numbers of generations !
                j_vec(2*D.num_stages) = j; 
                cur_ind = mysub2ind(D.num_params_vec, length(D.num_params_vec), j_vec);
                N_vec = demographic_parameters_to_n_vec(D, cur_ind); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation
                if(~filter_N_vec_internal(N_vec)) % get rid of this model 
                    max_j = j; 
                    break;
                end
            end
%             if(filter_N_vec_internal(N_vec)) % get rid of this model 
%                 continue;  % PROBLEM! this gets rid of a 'bad' model, and then will get rid of all sub-models, some of which might be 'good'
%             end
            [mu_vec_expansion_analytic] = FisherWright_Compute_SFS_Moments(N_vec, 0, num_moments); % compute moments with Formulas from Ewens
%            het_moment_mat_all_models(i,:) = mu_vec_expansion_analytic(:,end);
            
            for j=1:max_j % length(D.generations_vec{end}) % loop on all smaller numbers of generations !
                j_vec(2*D.num_stages) = j; 
                cur_ind = mysub2ind(D.num_params_vec, length(D.num_params_vec), j_vec);
                if(cur_ind == 250)
                    xxx = 24354325
                end
                het_moment_mat_all_models(cur_ind,:) = mu_vec_expansion_analytic(:,D.generations_vec{3}(j)+D.generations_vec{1}+D.generations_vec{2});
                init_N(cur_ind) = N_vec(1); % take first value of N 
            end            
        end
    end
    moments_time=cputime - moments_time
end


% need to convert sample to population!
[~, p_vec] = unique_with_counts(k_vec ./ n_vec); % get an empirical distribution 
total_sum = sum(p_vec);
p_vec = p_vec ./ sum(p_vec); % normalize

het_moment_mat_data = sum ( (k_vec./n_vec) .* (1-k_vec./n_vec) ) * L_correction_factor;  % compute empirical heterozygosity (corrected)

% het_moment_mat_data = moment_hist(x_vec, x_vec .* (1-x_vec) .* p_vec, 0, 0, 0) .* L_correction_factor;


epsilon = 0.05; % allow error in moments
region_het_moment_mat_all_models = het_moment_mat_all_models(:,1) .* 4 .* init_N .* mu; 
good_inds = find( (abs(region_het_moment_mat_all_models - het_moment_mat_data(:,1)) ./ het_moment_mat_data(:,1))  < epsilon );


i_ctr=1; log_like_mat = zeros(D.num_params, 1); compute_time=zeros(D.num_params, 1); 
for i=vec2row(good_inds) % 1:D.num_params
    
    length(D) % Loop on D: enumerate on many different demographic parameters    
    D.index = i;
    N_vec = demographic_parameters_to_n_vec(D, D.index); % D.generations, D.expan_rate, D.init_pop_size); % compute population size at each generation

    if(filter_N_vec_internal(N_vec))
        continue; 
    end
    
%    [demographic_x_vec, demographic_f_vec, ~, ~, demographic_compute_time] = ...
%        compute_allele_freq_spectrum_from_demographic_model(D, 0, compute_flag); % Try a grid of different values
    
    % Compute likelihood. This is trivial one-class likelihood (no mixture bullshit) so should be fast !!!
    rare_cumulative_per_gene = []; % set dummy variables
    target_size_by_class_vec = [mu, mu, mu]; % [neutral, null, missense]
    full_flag = 0; % use summary statistics 
    null_w_vec = 1; % NULL_C % assume all alleles are 'null' (but s=0 so actually neutral) 
    [log_like_mat(i), ~, compute_time(i)] = ... % compute likelihood (here vary only alpha)
        compute_two_class_log_likelihood(s, alpha, beta, rare_cumulative_per_gene, target_size_by_class_vec, D, ...
        X, [], [], null_w_vec, ...
        0, full_flag, []); % don't include phenotype !!
    fprintf('run good ind %ld out of %ld, ', i_ctr, length(good_inds));
    fprintf(' Loglike=%f, cur-time=%f, total-time=%f\n', log_like_mat(i), compute_time(i), sum(compute_time(1:i)));
    
    i_ctr=i_ctr+1;
end

[max_LL, max_J] = max(log_like_mat(good_inds)); max_J = good_inds(max_J); % maximize likelihood. (Take only good inds)
D = D{max_J};


% New: allow D_opt to use the fastneutrino method of Song et al.: 
switch D_opt 
    case 'neutrino'
        % prepare output file 
        
        fit_run_str = 'fastNeutrino'; % ...; % generate running command 
        
end

% Internal function for filtering out obviously unrealistic demographic models 
function filter_flag = filter_N_vec_internal(N_vec)    % Filter first unreasonable models !!!
filter_flag=0;
if(max(N_vec) > 10^9) % 1 billion
    filter_flag=1;
end
if(length(N_vec) > 4000) % we model maximum of 4000 generations
    filter_flag=1;
end


