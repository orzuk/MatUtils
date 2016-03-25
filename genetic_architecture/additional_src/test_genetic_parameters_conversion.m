% Convert between different quantaties and make sure conversion works
% function test_genetic_parameters_conversion()
AssignGeneralConstants;
test_back_and_forth = 1; test_lambda_s = 0;
h_scale = 'binary'
iters = 10; % draw different parameters
f_vec = rand(iters, 1) ./ 2; % assume always RAF < 0.5
mu = rand(iters, 1) ./ 5; % assume prevalence < 20%
h_add_input = rand(iters, 1) ./ 3; % assume heritability up to 33% %  h_add_input = 0.03;



mu=0.01; f = 0.001; 
beta_vec =  0:0.005:20; %   logspace(-5, 1, 300)
h_liab_vec = 2.*f.*(1-f).*beta_vec.^2; 

for i=1:length(beta_vec)
%     [grr_vec(i), grr_vec_binary(i)] = beta_to_genetic_relative_risk(beta_vec(i), f, mu, 'binary'); 
%     [~, ~, ~, h_disease_vec_binary(i)] = genetic_relative_risk_to_heritability(f, grr_vec_binary(i), mu, 'binary');
%     h_disease_vec_binary(i) = 2*h_disease_vec_binary(i); 
%     [~, ~, ~, h_disease_vec(i)] = genetic_relative_risk_to_heritability(f, grr_vec(i), mu, 'binary');
%     h_disease_vec(i) = 2*h_disease_vec(i); 
    h_disease_vec(i) = heritability_scale_change(h_liab_vec(i), 'binary', mu); 
    [h_disease_vec_binary(i) rr_vec(i)] = beta_to_variance_explained(beta_vec(i), f, [], 'binary', [], 'binary', mu); 
    h_disease_vec_binary(i) = 2 .* h_disease_vec_binary(i); 
end

bad_inds = find(h_liab_vec > 1);  % bad_inds = find(grr_vec > 1/mu); 
good_inds = setdiff(1:length(beta_vec), bad_inds); % Don't allow grr to be higher than 1 / mu


LT_dir = '../../common_disease_model/docs/figs/LT_big_effects/';    
figure; hold on; plot(h_liab_vec(good_inds), h_disease_vec(good_inds)); 
plot(h_liab_vec(good_inds), h_disease_vec_binary(good_inds), 'r'); xlim([0 0.1]);
xlabel('Var. on liab. scale');
ylabel('Var. on disease scale'); 
legend({'naive-LT', 'MoG-transofrm'}); legend('boxoff');
my_saveas(gcf, fullfile(LT_dir, 'var_explained_disease_vs_liability'), {'epsc', 'pdf'}); 


figure; hold on;  plot(beta_vec(good_inds), h_disease_vec(good_inds) ./ h_liab_vec(good_inds)); 
plot(beta_vec(good_inds), h_disease_vec_binary(good_inds) ./ h_liab_vec(good_inds), 'r'); 
legend({'naive-LT', 'MoG-transofrm'}); legend('boxoff');
xlabel('\beta'); ylabel('Inflation factor of var. explained  h_{disease}^2 / h_{liab}^2 '); % xlim([0 2]);
my_saveas(gcf, fullfile(LT_dir, 'beta_vs_var_explained_inflation_factor'), {'epsc', 'pdf'}); 


figure; hold on; plot(rr_vec(good_inds), h_disease_vec_binary(good_inds)); 
xlabel('RR'); ylabel('Var. on disease scale'); 
my_saveas(gcf, fullfile(LT_dir, 'rr_vs_var_explained_disease'), {'epsc', 'pdf'}); 

figure; hold on; plot(rr_vec(good_inds),h_disease_vec_binary(good_inds) ./ h_liab_vec(good_inds)); 
xlabel('RR'); ylabel('Inflation factor of var. explained h_{disease}^2 / h_{liab}^2 ');
my_saveas(gcf, fullfile(LT_dir, 'rr_vs_var_explained_inflation_factor'), {'epsc', 'pdf'}); 

figure; semilogx(rr_vec(good_inds)-1,h_disease_vec_binary(good_inds) ./ h_liab_vec(good_inds));
xlabel('RR-1 (log-scale)'); ylabel('Inflation factor of var. explained h_{disease}^2 / h_{liab}^2 ');
my_saveas(gcf, fullfile(LT_dir, 'rr_vs_var_explained_inflation_factor_log'), {'epsc', 'pdf'}); 


if(test_back_and_forth) % Here convert between genetic relative risk and heritability 
    for i=1:iters
        GRR(i) = heritability_to_genetic_relative_risk(h_add_input(i), h_scale, f_vec(i), mu(i));
        h_liab_input(i) = heritability_scale_change(h_add_input(i), 'liability', mu(i));
        beta_input(i) = sqrt(h_liab_input(i) ./ (f_vec(i) .* (1-f_vec(i)))); 

        
        GRR(i) = heritability_to_genetic_relative_risk(h_add_input(i), h_scale, f_vec(i), mu(i));
        grr_from_beta(i) = beta_to_genetic_relative_risk(beta_input(i), f_vec(i), mu(i), 'binary'); 
        beta_from_grr(i) = genetic_relative_risk_to_beta(f_vec(i), grr_from_beta(i), mu(i), 'binary'); 
        
        [lambda_s_vec lambda_s_add(i) lambda_mz_add(i) h_add(i) V_add ...
            lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab(i)] = ...
            genetic_relative_risk_to_heritability(f_vec(i), GRR(i), mu(i));
        
        h_liab2(i) = genetic_relative_risk_to_variance_explained(f_vec(i), GRR(i), mu(i), 'binary')
        %    genetic_relative_risk_to_heritability([f_vec f_vec]', [GRR GRR]', mu)
        should_be_zero = h_add(i) - h_add_input(i)
        should_be_zero_liability_scale = h_liab(i) - h_liab_input(i)
        should_be_zero_beta = beta_from_grr(i) - beta_input(i) 
    end
    
    parameters_labels = {'RAF', 'Prevalence', 'GRR', 'h^2-binary', 'h^2-liability', 'lambda_{MZ}', 'lambda_{s}'};
    parameters_data = [f_vec mu vec2column(GRR) h_add' h_liab' lambda_s_add' lambda_mz_add'];
    parameters_data = str2num(num2str(parameters_data, 4)); % round
    
    cell_iter_vec = num2cell(1:iters)'; % save inputs&outputs in a unit-tests file
    save_expression_mat_as_text_generic(parameters_labels, cell_iter_vec, parameters_data, ...
        '../../common_disease_model/data/test_genetic_parameters_conversion.txt', 1);
    
    
    % Eliana's example:
    [lambda_s_vec ...
        lambda_s_add lambda_mz_add h_add V_add ...
        lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab] = ...
        genetic_relative_risk_to_heritability(0.2, 2, 0.05)
    should_be_GRR_2 = heritability_to_genetic_relative_risk(h_liab, 'liability', 0.2, 0.05)
    %heritability.to.grr(grr.to.varexpl(0.2,2,0.05),"liability",0.2,0.05)
    
    GRR = 1.2; f = 0.2; mu = 1/9;
    [lambda_s_vec ...
        lambda_s_add lambda_mz_add h_add V_add ...
        lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab] = ...
        genetic_relative_risk_to_heritability(f, GRR, mu)
    should_be_GRR_1_2 = heritability_to_genetic_relative_risk(h_liab, 'liability', f, mu)
end % test back and forth


if(test_lambda_s) % Test conversion of liability to lambda_z:
    tol = 0.05; h_liab_vec = tol:tol:1; mu = 0.01; f = 0.5;
    
    N = 50; % 10; % large number of loci -> good Gaussian approximation
    compute_method_flag = 'sampling'; % 'enumerate'; % 'sampling'; % 'enumerate'; % 'sampling'
    iters = 1000; arch_str = 'circuit'; n_gates = N;
    max_generations = 1; f_vec = repmat(f, N, 1);  % N loci
    num_gates = N+1; % N+5;
    circuit = zeros(n_gates);
    % circuit(1,N+1) = 1;  circuit(2,N+1) = 1;
    for i=1:N
        circuit(i,N+1) = 1;
        %    circuit(i+2, N+i+1) = 1;
    end
    
        
    circuit(N+1,N+1) = LIABILITY; % just add loci contribution
    
    figure; graph_draw(circuit);
    circuit_gate_params = cell(num_gates,1);
    
    params_struct = [];
    params_struct.z_std = 0; % the environmental st.d.
    params_struct.min_freq = 0;
    params_struct.max_freq = 1;
    params_struct.N = N;
    %    params_struct.linear_coef_vec = log(GRR_marginal);
    params_struct.circuit = circuit;
    
    mu_vec = [0.1 0.05 0.1 0.2];
    for j=1:length(mu_vec)
        mu = mu_vec(j);
        threshold = norminv(1 - mu);
        for i=5:length(h_liab_vec) % here create liability architecture and compute empirically
            x = rand(N, iters) < f; f_vec = repmat(f, 1, N); % randomize each time (avoid biases in x)
            %    params_struct = set_architecture_params
            h_liab = h_liab_vec(i);
            circuit_gate_params{N+1}.a = h_liab;
            circuit_gate_params{N+1}.b = threshold + sqrt(h_liab)*sqrt(N) * sqrt(f/(1-f));
            circuit_gate_params{N+1}.linear_coef_vec = repmat(sqrt(h_liab) / sqrt(N*f*(1-f)), 1, N);
            
            %         circuit_gate_params{N+2}.a = h_liab*2 / sqrt(N); circuit_gate_params{N+2}.b = -h_liab*sqrt(N); % set affine transfomration
            %         circuit_gate_params{N+3}.a = sqrt(1-h_liab^2);
            %         circuit_gate_params{N+5}.b = threshold;
            params_struct.circuit_gate_params = circuit_gate_params;
            
            [mu_computed(i) V(i) v_environment(i) v_genetic(i) mz_twin_risk(i) H_empirical(i)] = ... % (can't compute liability yet) h_liability] = ...
                compute_architecture_statistics('circuit', ...
                f_vec, params_struct, [], [], iters, ... % one output per x
                compute_method_flag) % compute several moments and other stuff for architecture
            [v_marginal(:,i) GRR_marginal(:,i) p_z_x_marginal h_add(i) h_liability(i)] = ...
                compute_architecture_statistics_marginal('circuit', ...
                f_vec, params_struct, [], [], iters, ...
                compute_method_flag, mu_computed(i));
            
            
            relative_sampling_iters = 1000; % how many families to simulate
            [family_risk{i} family_tree relative_risk{i} sibs_genotype{i} sibs_freqs{i} ...
                H_from_twins(i) h_add_from_twins(i) h_liability_from_twins(i) h_liability_from_twins_ADE(i)] = ...
                compute_architecture_family_risk('circuit', ...
                f_vec, params_struct, relative_sampling_iters, ...
                compute_method_flag, 1);
            lambda_s{j}(i) = family_risk{i}(end-1) / mu_computed(i);  % compute lambda_s
            lambda_mz{j}(i) =  family_risk{i}(end)/ mu_computed(i)
            
            if(N == 10)
                [mu_computed_enumerate(i) V_enumerate v_environment_enumerate v_genetic_enumerate ...
                    mz_twin_risk_enumerate h_empirical_enumerate(i)] = ... % (can't compute liability yet) h_liability] = ...
                    compute_architecture_statistics('circuit', ...
                    f_vec, params_struct, [], [], iters, ... % one output per x
                    'enumerate') % compute several moments and other stuff for architecture
                [v_marginal_enumerate GRR_marginal_enumerate p_z_x_marginal_enumerate h_add_enumerate(i)] = ...
                    compute_architecture_statistics_marginal('circuit', ...
                    f_vec, params_struct, [], [], iters, ...
                    'enumerate', mu_computed_enumerate(i)); % Important to take the computed mu (not theoretical)
            end
            
            z = apply_circuit(x', circuit, circuit_gate_params);
            mu_vec_empirical(i) = mean(z);
        end
        figure;plot(h_liab_vec, mu_vec_empirical, '.'); hold on;
        plot(h_liab_vec, mu_computed, 'g.'); hold on;
        plot(h_liab_vec, repmat(mu, length(h_liab_vec),1), 'r');
        if(N == 10)
            plot(h_liab_vec, mu_computed_enumerate, 'c.'); hold on;
        end
        title(['empirical prevalence as function of heritability (should be constant) \mu=' num2str(mu)]);
        legend('sampled', 'sampled-computed', 'true'); %  'enumerate', );
        
        [h_vec h_broad_vec h_vec_exact] = heritability_scale_change(h_liab_vec, 'binary', mu);
        figure; hold on; plot(h_liab_vec, h_vec, '.');
        hold on; plot(h_liab_vec, h_vec_exact, 'r.');
        plot(h_liab_vec, h_broad_vec, 'gx');
        plot(h_liab_vec, h_add, 'm.');
        plot(h_liab_vec, h_empirical, 'g--');
        if(N == 10)
            plot(h_liab_vec, h_add_enumerate, 'kx');
            plot(h_liab_vec, h_empirical_enumerate, 'g+');
        end
        
        legend_vec = {'approx.', 'exact', 'broad', 'empiric', 'broad-empiric'};
        if(N == 10)
            legend_vec = [legend_vec {'empiric-enumerate', 'broad-empiric-enumerate'}]; % Problem: approx and exact are way off!!!!
        end
        legend(legend_vec);
        figure; plot(h_liab_vec, h_vec - h_vec_exact, '.'); title('approximation quality');
        
        
    end % loop on mu
    
    figure; hold on; legend_vec = {}; % plot lambda_s vs. heritability: empirical and theoretical  
    for j=1:length(mu_vec) % plot for different mu's
        mu = mu_vec(j);
        plot(h_liab_vec, log(lambda_s{j}), [color_vec(j) '.']);
        lambda_s_vec_theoretical{j} = heritability_to_sibling_relative_risk(h_liab_vec, 'liability', mu);
        plot(h_liab_vec, log(lambda_s_vec_theoretical{j}), color_vec(j));
        legend_vec{2*j-1} = ['\mu=' num2str(mu) '(emp.)']; legend_vec{2*j} = ['\mu=' num2str(mu) '(anal.)'];
        
    end
    ylabels = get(gca, 'YTickLabel');
    ylabels = num2str(exp(str2num(ylabels)), 3);
    set(gca, 'YTickLabel', ylabels);
    title(['Heritability (liability-scale) vs. \lambda_s']);
    xlabel('h_{liab}^2'); ylabel('\lambda_s');
    legend(legend_vec);

    figure; hold on; legend_vec = {}; %  Plot the same for lambda_MZ
    for j=1:length(mu_vec) % plot for different mu's
        mu = mu_vec(j);
        plot(h_liab_vec, log(lambda_mz{j}), [color_vec(j) '.']);
        [lambda_s_vec_theoretical{j} lambda_mz_vec_theoretical{j}] = ...
            heritability_to_sibling_relative_risk(h_liab_vec, 'liability', mu);
        plot(h_liab_vec, log(lambda_mz_vec_theoretical{j}), color_vec(j));
        legend_vec{2*j-1} = ['\mu=' num2str(mu) '(emp.)']; legend_vec{2*j} = ['\mu=' num2str(mu) '(anal.)'];
        
    end
    ylabels = get(gca, 'YTickLabel');
    ylabels = num2str(exp(str2num(ylabels)), 3);
    set(gca, 'YTickLabel', ylabels);
    title(['Heritability (liability-scale) vs. \lambda_{MZ}']);
    xlabel('h_{liab}^2'); ylabel('\lambda_{MZ}');
    legend(legend_vec);

end % test lambda_s 
