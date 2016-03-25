% Run IBD estimation
AssignGeneralConstants;
queue_str = 'priority';
memory_needed = [2 3 10]; % memory required (in GB. 4 should be enough for 5000.)

%num_people_vec = {'500:500:1000', 
num_people_vec = {'500:500:2000', '500:500:4000', '500:500:10000'};
LP_model_str = {'linear', 'debug', 'P*', 'P*_full'};
num_snps = 1000;

% Qatar example
mean_Q = 0.0358; std_Q = 0.0395;
num_founders_Q = round(1/mean_Q);
num_snps = 1000;
num_generations_Q = 23;
IBD_iters = 100;
num_founders_vec = [10 50 num_founders_Q];
num_generations_vec = num_generations_Q;


if(machine == UNIX) % launch jobs
    for i_f = 1:length(num_founders_vec) % loop on # founders (determines k0)
        for i_g = 1:length(num_generations_vec) % loop on # of generations
            for j=1:length(num_people_vec) % loop on # people to run
                for i=1:length(LP_model_str) % loop on model
                    run_str = ['LP_model_str=''' LP_model_str{i} '''; num_snps=' num2str(num_snps) ...
                        '; num_people_vec=' num_people_vec{j} '; num_founders=' num2str(num_founders_vec(i_f)) ...
                        '; test_estimate_IBD_sharing_from_snps;']
                    log_file = ['out/IBD_sharing_' strrep(LP_model_str{i}, '*', '_star') ...
                        '_num_snps_' num2str(num_snps) ...
                        '_num_people_' strrep(num_people_vec{j}, ':', '_') '_num_founders_' ...
                        num2str(num_founders_vec(i_f)) '.out']
                    SubmitMatlabJobToFarm(run_str, log_file, queue_str,[], [], memory_needed(j)); % ,[], memory_needed);
                end
            end
        end
    end
end

% LP 10 model doesn't work !!!! 
LP_model_str = 'P*'; % 'debug'; % 'P*'; % 'debug'; % 'P*';
num_people_vec = [400:400:800]; % 1200]; % change number of people 
num_founders = num_founders_Q; % change to include only qatar example
num_generations = num_generations_Q;
IBD_iters = 1000; % # of iterations 
trait_type = 'disease';  %'quantitative'; % 'disease'; % 'disease'; % 'quantitative'; % 'disease'; % 'disease'
%num_founders = 50; num_generations = 5; % take k0 close to zero but allow lots of variation  
test_estimate_IBD_sharing_from_snps



return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compute_eliana_example = 0;
if(compute_eliana_example)
    
    genetic_relative_risk_to_familial_risk(0.3, 1.1) * ...
        genetic_relative_risk_to_familial_risk(0.1, 1.3)
    
    max_generations = 1; f_vec = repmat(f, N, 1);  % N loci    
    num_gates = 3; % N+5;
    circuit = zeros(3);
    % circuit(1,N+1) = 1;  circuit(2,N+1) = 1;
    for i=1:2
        circuit(i,3) = 1;
        %    circuit(i+2, N+i+1) = 1;
    end
    circuit(3,3) = AND; % just add loci contribution
    
    figure; graph_draw(circuit);
    params_struct.circuit = circuit;
    params_struct.z_std = 0; % the environmental st.d.
    params_struct.min_freq = 0.1;
    params_struct.max_freq = 0.2;
    params_struct.N = 2;
    
    %params_struct = [];
    [lambda_vec family_tree relative_risk sibs_genotype sibs_freqs] = ...
        compute_architecture_family_risk('circuit', ...
        [0.1 0.3]', params_struct, [], ...
        'enumerate', 1);    
    
    mu = 0.1 * (1-0.3*0.1) + 0.2 * 0.3*0.1
    prob_2 = [0.9*0.05*0.7*0.15 + 0.9*0.05*0.3*0.65 + 0.1*0.55*0.7*0.15]*2;
    prob_4 = 0.1*0.55*0.3*0.65;
    
    prob_1 = 1-prob_2-prob_4    
    lambda_s = (prob_1 *0.01 + prob_2 * 0.02 + prob_4 * 0.04) / mu^2




    % Test with actual data
    for j=1:100
        iters=1000000;
        shared_genotype_mat = rand(iters,2)> 0.5; % set if we share each genotype
        genotype_mat = rand(iters,4) < repmat([0.1 0.3 0.1 0.3], iters,1);
        genotype_mat(:,3) = shared_genotype_mat(:,1) .* genotype_mat(:,1) + ...
            (1-shared_genotype_mat(:,1)) .* genotype_mat(:,3);
        genotype_mat(:,4) = shared_genotype_mat(:,2) .* genotype_mat(:,2) + ...
            (1-shared_genotype_mat(:,2)) .* genotype_mat(:,4);
        
        base_mu=0.4;
        risk_vec1 = sum(genotype_mat(:,1:2),2) == 2;
        risk_vec2 = sum(genotype_mat(:,3:4),2) == 2;
        
        z_vec1 = rand(iters,1) < base_mu .* (1+risk_vec1);
        z_vec2 = rand(iters,1) < base_mu .* (1+risk_vec2);
        mean(z_vec1)
        mean(z_vec2)
        
        lambda_s(j) = mean(z_vec1 .* z_vec2) ./ (mean(z_vec1)*mean(z_vec2))
    end % loop to compute sib-risk 

end % compute eliana example 



try_heritability_conversion=1;
if(try_heritability_conversion)
    
    % Example: One big locus. Additive:
    n=1000000;
    beta = sqrt(2); 
    g = beta .* ((rand(n, 1) > 0.5)-0.5); 
    gg = randn(n, 1) .* sqrt(0.5);
    e = randn(n, 1) .* sqrt(0.5);
    x = g+e; % liability
    xx = gg+e;
    mu = 0.1; x_mu = norminv(1-mu);
    z = x>x_mu; % disease 
    zz = xx>x_mu;
    h_liab = corr(g, x)^2
    h_01 = corr(g, z)^2
    heritability_scale_change(h_01, 'liability', mu)

    h_liab_xx = corr(gg, xx)^2
    h_01_xx = corr(gg, zz)^2
    heritability_scale_change(h_01_xx, 'liability', mu)

        
    
    
    % Try log-normal distribution: small effect. Non-additive
    beta = 0.1; 
    g = beta .* ((rand(n, 1) > 0.5)-0.5); 
    e = randn(n, 1) .* sqrt(1-h_liab);
    x_gaussian = g+e; mean(x_gaussian), std(x_gaussian)
    x = exp(g+e); % liability
    x_mu = exp(norminv(1-mu))
    z = x>x_mu; % disease 
    prevlaence_is = mean(z)
    prevalence_plus_is = mean(z(g>0))
    prevalence_minus_is = mean(z(g<0))

    
    mean(x), should_be_mean = exp(0 + 1/2)
    var(x), should_be_var = (exp(1)-1)*exp(1)
    
    h_liab = corr(g, min(x,10))^2
    h_liab = corr(g, x)^2
    true_beta_liab = 2*sqrt(h_liab)
    h_01 = corr(g, z)^2 
    standard_Demptseter_should_be_h_liab = heritability_scale_change(h_01, 'liability', mu)
    standard_Demptseter_should_be_h_liab2 = h_01 * mu*(1-mu) / ...
        normpdf(norminv(1-mu))^2  

    
    
    heritability_scale_change(h_01, 'liability', mu)
    generalized_Dempster_should_be_h_liab = h_01 * mu*(1-mu) / ...
        (exp(-log(x_mu)^2/2)/(x_mu*sqrt(2*pi)))^2  
    ratio_is = generalized_Dempster_should_be_h_liab/ h_liab

    figure; hold on; hist_density(x, 1000); xlim([0 10]); 
    x_vec = 0:0.01:10; plot(x_vec, (exp(-log(x_vec).^2./2)./(x_vec.*sqrt(2*pi))), 'r'); 



    delta_prevalence = prevalence_plus_is-prevlaence_is
    delta_prevalence_predicted = true_beta_liab * 0.5 * exp(-log(x_mu)^2/2)/(x_mu*sqrt(2*pi))


end

