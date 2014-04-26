% Compute fraction of nulls given that allele frequency is below f (for equilibrium)
%
% Input:
% s - selection coefficient
% N - effective population size
% alpha_vec - fraction of null alleles at birth
% f_vec - top allele frequency (threshold)
% demographic_models_struct - model parameters and distribution for null
% and neutral alleles
% model_str - name of model
% n - sample size (we assume we look at a population. For case-control formula is different)
% compute_method - analytic (default) or simulations
%
% Output:
% rho_vec - fraction of null alleles below each allele frequency, when we
% use weighted distribution \psi_s(f) (weighted by allele frequency)
%
function rho_vec = frac_null_conditional_on_freq_less_f(s, N, alpha_vec, f_vec, ...
    demographic_models_struct, model_str, n, polyphen_error, compute_method, iters)

if(~exist('demographic_models_struct', 'var') || isempty(demographic_models_struct))
    demographic_models_struct = []; % set default
end
if(~exist('compute_method', 'var') || isempty(compute_method))
    compute_method = 'analytic'; % set default
end

if(~exist('iters', 'var') || isempty(iters))
    iters = 10000; % number of alleles
end

num_thresholds = length(f_vec); % get number of different thresholds 

if(isempty(demographic_models_struct))  % here compute analytically for equilibrium
    S=4.*N.*s; % get effective selection coefficient
    rho_vec = alpha_vec .* (phi_s_integral(f_vec, S, 1) - phi_s_integral(0, S, 1));
    rho_vec = rho_vec ./ (rho_vec + (1-alpha_vec) .* (phi_s_integral(f_vec, 0, 1) - phi_s_integral(0, 0, 1)));
else % here compute from empirical distribution
    model_ind = strcmp(model_str, demographic_models_struct.model_str);
    [~, s_ind] = min(abs(s - demographic_models_struct.data{model_ind}.s_vec)); % take closest s
    
    %%%    s_ind=2; % TEMP BAD!!
    
    %    s_ind = length(demographic_models_struct.data{model_ind}.s_vec)+1-s_ind; % take reverse (p_vec is in reversed order)
    psi_s_cumulative = demographic_models_struct.data{model_ind}.p_vec(s_ind,:);
    psi_neutral_cumulative = demographic_models_struct.data{model_ind}.p_vec(end,:);
    
    if(~exist('polyphen_error') || isempty(polyphen_error))
        gamma_null = 1; gamma_neutral = 1;
    else
        gamma_null = 1-polyphen_error(1); gamma_neutral = polyphen_error(2);
    end
    
    if(~exist('n', 'var') || isempty(n)) % compute rho for population frequency
        rho_vec = alpha_vec .* gamma_null .* psi_s_cumulative ./ ...
            (alpha_vec .* gamma_null .* psi_s_cumulative + (1-alpha_vec) .* gamma_neutral .* psi_neutral_cumulative);
    else  % compute rho for sample frequency. n cannot be too large
        
        switch compute_method
            case 'analytic'
                k_vec = round(f_vec .* 2*n); max_k = max(k_vec); % number of carriers
                rho_vec = zeros(length(k_vec), 1);
                x_vec = demographic_models_struct.data{model_ind}.x_vec; % allele-frequency values
                num_x = length(x_vec);
                binom_mat = zeros(num_x, max_k+1); % binom_mat(true-freq., sample-freq.) = Pr(sample-freq. | true-freq.)
                for i=1:(max_k+1)
                    binom_mat(:,i) = binopdf(i-1, 2*n, x_vec); % Pr(sample-freq=(i-1)/2n | true-freq = x_vec)
                end
                psi_s_density = [psi_s_cumulative(1) diff(psi_s_cumulative)];
                psi_neutral_density = [psi_neutral_cumulative(1) diff(psi_neutral_cumulative)]; % Get density at each frequency bin
                
                % % % OLD WAY (has a bug!!)
                % % %                 for i=1:length(k_vec) % loop on thresholds
                % % %                     if(k_vec(i) == 0) % threshold is lower than 1/2n. We don't 'pass' any allele
                % % %                         rho_vec(i) = alpha_vec;
                % % %                     else
                % % %                         binom_vec = vec2row(sum( (0:k_vec(i))' .* binom_mat(:, 1:k_vec(i)+1), 2)); % sum over all TRUE frequencies giving the same sample frequency. Why not from one?
                % % % %                                                  rho_vec(i) = sum(alpha_vec .* gamma_null .* psi_s_density .* binom_vec) ./ ... % sum over all true frequencies
                % % % %                                                      sum((alpha_vec .* gamma_null .* psi_s_density + ...
                % % % %                                                      (1-alpha_vec) .* gamma_neutral .* psi_neutral_density) .* binom_vec);
                % % %                          rho_vec(i) = sum(alpha_vec .* gamma_null .* psi_s_density .* binom_vec .* (i./x_vec')) ./ ... % sum over all true frequencies
                % % %                              sum((alpha_vec .* gamma_null .* psi_s_density   + ...
                % % %                              (1-alpha_vec) .* gamma_neutral .* psi_neutral_density) .* binom_vec .* (i./x_vec'));
                % % %                     end
                % % %                 end
                
                for i=1:(max_k+1) % loop on SAMPLE frequencies
                    rho_s_sample_density_vec(i) =  sum(  (psi_s_density ./ x_vec') .* binom_mat(:,i)' ); % New: sum the other way (over true frequencies
                    rho_neutral_sample_density_vec(i) = sum(  (psi_neutral_density ./ x_vec') .* binom_mat(:,i)' );
                end
                
                rho_s_sample_density_vec = cumsum(rho_s_sample_density_vec .* ((0:max_k) ./ (2*n))); % multiply by alleles
                rho_neutral_sample_density_vec = cumsum(rho_neutral_sample_density_vec .* ((0:max_k) ./ (2*n))); % multiply by alleles
                
                for i=1:length(k_vec) % loop on k % k_vec
                    run_i = i
                    I_sample = find((0:2*n)./(2*n) <= f_vec(i), 1, 'last'); % IMPORTANT! TAKE INDEX BEFORE OR AFTER!
%                    I_sample = min(I_sample, max(k_vec)+1);
                    rho_s_density_vec(i) = rho_s_sample_density_vec(I_sample);
                    rho_neutral_density_vec(i) = rho_neutral_sample_density_vec(I_sample);
                end
                
                rho_vec = alpha_vec .* gamma_null .* rho_s_density_vec ./ ...
                    ( alpha_vec .* gamma_null .* rho_s_density_vec + ...
                    (1-alpha_vec) .* gamma_neutral .* rho_neutral_density_vec ); % set rho
            case 'sampling' % NEW! Here compute by sampling alleles! (should give the same results)
                [class_vec f_pop_vec f_sample_vec] = ...
                    simulate_allele_frequencies(iters, demographic_models_struct, alpha_vec, n);
                
                x_vec = demographic_models_struct.data{1}.x_vec;
                for i=1:length(x_vec)
                    run_i = i
                    total_null_freq(i) = sum(f_sample_vec(f_sample_vec < x_vec(i)) .* class_vec(f_sample_vec < x_vec(i)));
                    total_neutral_freq(i) = sum(f_sample_vec(f_sample_vec < x_vec(i)) .* (1-class_vec(f_sample_vec < x_vec(i))));
                    
                    total_null_freq_pop(i) = sum(f_pop_vec(f_pop_vec < x_vec(i)) .* class_vec(f_pop_vec < x_vec(i)));
                    total_neutral_freq_pop(i) = sum(f_pop_vec(f_pop_vec < x_vec(i)) .* (1-class_vec(f_pop_vec < x_vec(i))));
                end
                %                 total_null_freq = total_null_freq .* demographic_models_struct.data{1}.p_vec(1,end);
                %                 total_neutral_freq = total_neutral_freq .* demographic_models_struct.data{1}.p_vec(2,end);
                rho_vec = total_null_freq ./ (total_null_freq + total_neutral_freq);
                rho_vec_pop = total_null_freq_pop ./ (total_null_freq_pop + total_neutral_freq_pop);
                
                %
        end % switch compute method
        
    end % if finite sample (n exists)
end % if equilibrium (or demographic models from simulations)


% Internal function: simulate allele frequencies
%

% Output:
% f_vec - vector of allele frequencies
% class_vec - vector of binary variables saying if allele is NULL (1) or NEUTRAL (0)
function [class_vec f_pop_vec f_sample_vec] = ...
    simulate_allele_frequencies(iters, demographic_models_struct, alpha, n)


class_vec = zeros(iters,1); f_pop_vec = zeros(iters,1);

num_null = binornd(iters, alpha); num_neutral = iters-num_null;
class_vec(1:num_null) = 1; % First simulate class

x_vec = demographic_models_struct.data{1}.x_vec; % Simulate from population (population frequencies)

p_vec_density = [0 diff(demographic_models_struct.data{1}.p_vec(1,:))];
p_vec_density = p_vec_density ./ x_vec'; % NEW! get-rid of frequency weighting
w_null = weighted_rand(p_vec_density, num_null);  % Next simulate allele frequencies
f_pop_vec(1:num_null) = x_vec(w_null);
p_vec_density = [0 diff(demographic_models_struct.data{1}.p_vec(2,:))];
p_vec_density = p_vec_density ./ x_vec'; % NEW! get-rid of frequency weighting
w_neutral = weighted_rand(p_vec_density, num_neutral);
f_pop_vec((1+num_null):iters) = x_vec(w_neutral);

f_sample_vec = binornd(2*n, f_pop_vec) ./ (2*n); % randomize sample frequencies

