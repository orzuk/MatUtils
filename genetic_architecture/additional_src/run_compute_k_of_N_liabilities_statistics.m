AssignGeneralConstants();

% Produce some figures for k-of-N liabilities model
epistasis_figs_dir = '../../common_disease_model/figures/k_of_N_model';
arch_type = 'liability';  % 'pathways'; % % 'pathways';
S = {};
mu = 0.0001; max_N=10; % set maximal number of loci and prevalence
h_x = 1; % 0.8; % set heritability of each liability
%save_S = S; S= T;
format_fig_vec = 'jpg'

switch arch_type
    case 'liability'
        for N=1:max_N
            for k=1:1% N
                run_K_of_N = [k N]
                % Determine 'best' mu_l and h_x
                
                if(k == 1)
                    mu_l = 1 - (1-mu)^(1/N);
                else % this isn't very accurate ..
                    mu_l = fminbnd(@(x) abs(binocdf(k-1, N, x)-(1-mu)), 0, 1); % find mu_l that keeps the prevalence
                end
%                mu_l = 0.0001; 
                [~, S{N,k}] = compute_k_of_N_liabilities_statistics( N, k, mu_l, h_x); %
                S{N,k}.mu_l = mu_l;
            end
        end
        
        h_liab_vec = zeros(N);
        h_liab_from_twins_vec = zeros(N);
        lambda_s_vec = zeros(N);
        lambda_mz_vec = zeros(N);
        H_vec = zeros(N);
        h_vec = zeros(N);
        mu_vec = zeros(N);
        mu_l_vec = zeros(N);
        
        for N=1:max_N % collect results and plot figures
            for k=1:1% N
                h_liab_vec(N,k) = S{N,k}.h_liab;
                h_liab_from_twins_vec(N,k) = S{N,k}.h_liab_twins;
                H_vec(N,k) = S{N,k}.H;
                h_vec(N,k) = S{N,k}.h;
                lambda_s_vec(N,k) = S{N,k}.lambda_s;
                lambda_mz_vec(N,k) = S{N,k}.lambda_MZ;
                mu_vec(N,k) = S{N,k}.mu;
                mu_l_vec(N,k) = S{N,k}.mu_l;
            end
        end
        for N=1:max_N % collect results and plot figures
            figure; hold on;
            plot(1:N, h_liab_vec(N,1:N), '*'); plot(1:N, h_liab_vec(N,1:N));
            plot(1:N, h_liab_from_twins_vec(N,1:N), 'r*');
            plot(1:N, h_liab_from_twins_vec(N,1:N), 'r');
            title(['additive heritability (liability) vs. heritability estimated from twins (ACE). N=' ...
                num2str(N) ' \mu=' num2str(mu) ' h_x^2 = ' num2str(100*h_x) '%']);
            xlabel('K'); ylabel('h^2');
            legend({'h^2 fitted', '', 'h^2 from twin-study'}); set(gca, 'XTICK', 1:N)
            my_saveas(gcf, fullfile(epistasis_figs_dir, ...
                ['heritability_additive_estimated_vs_real_N_' num2str(N) '_mu_' num2str(mu) '_h_' num2str(h_x)]), format_fig_vec); % save figs.
            
            figure; hold on;
            plot(1:N, lambda_s_vec(N,1:N), '*');     plot(1:N, lambda_s_vec(N,1:N));
            plot(1:N, lambda_s_vec(N,1:N).^2, 'o');     plot(1:N, lambda_s_vec(N,1:N).^2, '--');
            plot(1:N, lambda_mz_vec(N,1:N), 'r*');   plot(1:N, lambda_mz_vec(N,1:N), 'r');
            title(['Sibling and mz twin risk. N=' num2str(N) ' \mu=' num2str(mu) ...
                ' h_x^2 = ' num2str(h_x) '%']);
            legend({'\lambda_s', '', '\lambda_s^2', '', '\lambda_{MZ}'}); set(gca, 'XTICK', 1:N)
            xlabel('K'); ylabel('familial risk');
            my_saveas(gcf, fullfile(epistasis_figs_dir, ...
                ['familial_risk_N_' num2str(N) '+mu_' num2str(mu) '_h_' num2str(h_x)]), format_fig_vec); % save figs.
            
            figure; hold on;
            plot(1:N, h_vec(N,1:N), '*'); plot(1:N, h_vec(N,1:N));
            plot(1:N, H_vec(N,1:N), 'r*');     plot(1:N, H_vec(N,1:N), 'r');
            title(['Broad vs. Narrow (additive) heritability on diseae (binary) scale. N=' ...
                num2str(N) ' \mu=' num2str(mu) ' h_x^2 = ' num2str(h_x) '%']);
            xlabel('K'); ylabel('H^2');
            legend({'h^2', '', 'H^2'}); set(gca, 'XTICK', 1:N)
            my_saveas(gcf, fullfile(epistasis_figs_dir, ...
                ['heritability_broad_vs_narrow_N_' num2str(N) '_mu_' ...
                num2str(mu) '_h_' num2str(h_x)]), format_fig_vec); % save figs.
            
            % %     figure; hold on;
            % %     plot(1:N, mu_vec, '*');
            % %     title('Prevalence of different models');
            % %     my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            % %         ['prevalence_N_' num2str(N)]), format_fig_vec); % save figs.
        end
        h_liab_overestimation_frac = h_liab_from_twins_vec ./ h_liab_vec;
        figure; imagesc_with_labels(h_liab_overestimation_frac, ...
            num2str_cell(num2cell(1:N)), num2str_cell(num2cell(1:N)),N+1); colorbar;
        xlabel('K'); ylabel('N');
        title(['Fraction over-estimated due to epistasis for different (K,N), \mu=' ...
            num2str(mu) ' h_x^2=' num2str(h_x)]);
        my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            ['epistasis_fraction_K_of_N_model']), format_fig_vec);
        h_liab_vec(h_liab_vec == 0) = NaN;
        figure; imagesc_with_labels(h_liab_vec, ...
            num2str_cell(num2cell(1:N)), num2str_cell(num2cell(1:N)),N+1); colorbar;
        xlabel('K'); ylabel('N');
        title(['Heritability of disease `(liab. scale) for different (K,N), \mu=' ...
            num2str(mu) ' h_x^2=' num2str(h_x)]);
        my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            ['disease_heritability_K_of_N_model']), format_fig_vec);
        
        
        res = 0.1;
        mu_vec = 10.^[-1:-1:-15]; %     res:res:0.5; % Here fix the architecture but change heritability and sib-risk
        h_x_vec = res:res:1-res; % fix heritability of each pathway
        N=2; K=1; MM = length(mu_vec); NN = length(h_x_vec);
        h_liab_mat = zeros(MM, NN);
        h_liab_from_twins_mat = zeros(MM, NN);
        H_mat = zeros(MM, NN);
        h_mat = zeros(MM, NN);
        lambda_s_mat = zeros(MM, NN);
        lambda_mz_mat = zeros(MM, NN);
        mu_mat = zeros(MM, NN);
        for i=1:length(mu_vec)
            do_i = i
            if(K == 1)
                mu_l_vec(i) = 1 - (1-mu_vec(i))^(1/N);
            else
                mu_l_vec(i) = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu_vec(i))), 0, 1); % find mu_l that keeps the prevalence
            end
            mu_l = mu_l_vec(i); 
            for j=1:length(h_x_vec)
                do_j = j
                [~, R] = compute_k_of_N_liabilities_statistics( N, K, mu_l, h_x_vec(j)); %
                h_liab_mat(i,j) = R.h_liab;
                h_liab_from_twins_mat(i,j) = R.h_liab_twins;
                H_mat(i,j) = R.H;
                h_mat(i,j) = R.h;
                lambda_s_mat(i,j) = R.lambda_s;
                lambda_mz_mat(i,j) = R.lambda_MZ;
                mu_mat(i,j) = R.mu;
            end
        end
        H_liab_epistasis = (h_liab_from_twins_mat - h_liab_mat) ./ h_liab_mat; % fraction of 'added variance'
        figure; imagesc_with_labels( H_liab_epistasis, ...
            num2str_cell(num2cell(h_x_vec)), num2str_cell(num2cell(mu_vec)), MM+1 ); colorbar;
        x_ticks = get(gca, 'XTickLabels')
        xlabel('one-liab-heritability (h_x^2)'); ylabel('prevalence (\mu)');
        title('Fraction over-estimated due to epistasis for different prevalence and genetic components');
        my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            ['epistasis_fraction_' num2str(K) '_of_' num2str(N) ...
            '_model_varying_prevalence_and_heritability']), format_fig_vec);
        %            ['heritability_additive_estimated_vs_real_N_' num2str(N) '_mu_' num2str(mu) '_h_' num2str(h_x)]), format_fig_vec); % save figs.
        
        
        figure; imagesc_with_labels( h_liab_mat, ...
            num2str_cell(num2cell(h_x_vec)), num2str_cell(num2cell(mu_vec)),  MM+1 ); colorbar;
        xlabel('one-liab-heritability (h_x^2)'); ylabel('prevalence (\mu)');
        title('Additive heritability assuming one liability scale for different prevalence and genetic components');
        my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            ['additive_heritability_assuming_one_liability_' num2str(K) '_of_' num2str(N) ...
            '_model_varying_prevalence_and_heritability']), format_fig_vec);
        
        figure; imagesc_with_labels(lambda_s_mat.^2 ./ lambda_mz_mat , ...
            num2str_cell(num2cell(h_x_vec)), num2str_cell(num2cell(mu_vec)),  MM+1 ); colorbar;
        xlabel('one-liab-heritability (h_x^2)'); ylabel('prevalence (\mu)');
        title('\lambda_s^2 / \lambda_{MZ} for different prevalence and genetic components');
        my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            ['lambda_s_vs_lambda_mz_' num2str(K) '_of_' num2str(N) ...
            '_model_varying_prevalence_and_heritability']), format_fig_vec);
        
        
        
        N=3; K=1; % Plot GRR for one vs. 3 liabilities
        freq_vec = [0.01 0.1 0.2]; mu_vec = [0.01 0.05 0.1];
        grr_vec = 1:0.01:3;
        symbol_vec = '.:-xo*';
        
        full_figure;  ctr=0; legend_vec = {}; h_liab_vec = {};
        for j=1:length(freq_vec)
            freq = freq_vec(j); % freq = 0.1;
            for mu_j = 1:length(mu_vec)
                ctr=ctr+1;
                mu_l = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu_vec(mu_j))), 0, 1);
                mu = sum(binopdf(K:N, N, mu_l));  %disease prevalence
                p_y_x_marginal = genetic_relative_risk_to_p_z_x_marginal(repmat(mu_l, 1, length(grr_vec)), ...
                    grr_vec', freq); % compute prob.
                p_total_tab = compute_k_of_N_joint_tab(N, K, mu_l); % Pr. one liability and three liabilities on)
                
                for i=1:length(p_y_x_marginal)
                    do_i = i
                    p_pathway_tab = vec2mat(p_y_x_marginal(i,:), 2)'
                    p_x_z(2,2) = p_pathway_tab(1,2)*p_total_tab(2,1) / sum(p_total_tab(:,1)) + ...
                        p_pathway_tab(2,2)*p_total_tab(2,2) / sum(p_total_tab(:,2));
                    p_x_z(2,1) = freq - p_x_z(2,2); % locus is on
                    p_x_z(1,2) = mu - p_x_z(2,2);
                    p_x_z(1,1) = 1 - p_x_z(2,2) - p_x_z(1,2) - p_x_z(2,1);
                    [locus_freq, disease_grr_vec(i), disease_mu] = ...
                        p_z_x_marginal_to_genetic_relative_risk(vec2row(mat2vec(p_x_z')));
                end
                plot(grr_vec, disease_grr_vec, [color_vec(j) symbol_vec(mu_j)]);
                legend_vec{ctr} = ['\mu = ' num2str(100*mu,3) '% f = ' num2str(freq,4)];
                
                [lambda_s_vec ...
                    lambda_s_add lambda_mz_add h_add V_add ...
                    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab h_liab_vec{ctr}] = ...
                    genetic_relative_risk_to_heritability(freq, grr_vec, mu_l) % Now plot # loci
                
            end % loop on prevalence
        end % loop on frequency
        xlabel('GRR for one liability'); ylabel('GRR for disease');
        title(['GRR for ' num2str(K) '-of-' num2str(N) ' liabilities model']);
        legend(legend_vec,1);
        axis([1 max(grr_vec) 1 max(grr_vec)]);
        my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            ['GRR_one_liability_vs_disease_' num2str(K) '_of_' num2str(N) ...
            '_model_varying_prevalence_and_heritability']), format_fig_vec);
        plot(grr_vec, 1 + (grr_vec-1)/3, 'k', 'linewidth', 3);
        full_figure;
        for j=1:ctr            % divide by two to move from alleles to loci
            plot(grr_vec, log10(1./(h_liab_vec{j}/2)),  [color_vec(j) symbol_vec(ceil(j/3))]);
        end
        xlabel('GRR for one liability'); ylabel('num loci for liability (log_{10})');
        title('GRR for vs. num. loci needed for liability model');
        legend(legend_vec,1);
        axis([1 max(grr_vec) 1 1/min(min_cell(h_liab_vec))]);
        my_saveas(gcf, fullfile(epistasis_figs_dir, ...
            'GRR_one_liability_vs_num_loci_needed_model_varying_prevalence_and_heritability'), ...
            format_fig_vec);
        
        
        
    case 'pathways'% Compute two-level model (pathway)
        max_M = 5; % num pathways
        max_L = 5; % num genes in a pathway
        freq_vec = [0.01 0.05 0.1 0.2]; % RAF
        grr = 2; % genetic relative risk for each locus
        mu = 0.01;
        %lambda_s = genetic_relative_risk_to_heritability(
        [lambda_s lambda_s_add lambda_mz XX YY lambda_s_mult lambda_mz_mult] = ...
            genetic_relative_risk_to_heritability([freq freq]', [grr grr]', mu)
        [lambda_s2 lambda_mz2] = ...
            genetic_relative_risk_to_familial_risk(freq, grr); % We don't really know the GRR for each locus !!!
        
        for f_ind = 1:length(freq_vec)
            freq = freq_vec(f_ind);
            init_mu = freq;
            [init_lambda_s init_lambda_mz] = genetic_relative_risk_to_familial_risk(freq, 1000000000000000000000);
            init_lambda_s = sqrt(init_lambda_s)
            init_lambda_mz = sqrt(init_lambda_mz)
            h_x = 1;
            for M=5:max_M % num pathways
                for L=5:max_L % num genes in a pathway
                    num_pathways_and_size = [M L]
                    N = M*L; % total number of loci
                    for K1 = 2:2 %2*L % Compute for pathway (need at least K1 to make pathway damaged)
                        lambda_s_pathway = architecture_k_of_N( 2*L, K1, freq, init_lambda_s);
                        lambda_mz_pathway = architecture_k_of_N( 2*L, K1, freq, init_lambda_mz);
                        pathway_freq = 1-binocdf(K1-1, 2*L, freq); % probability a pathway is 'affected'
                        for K2 = 5:M % need at least K2 to pathways to make disease risk high
                            lambda_s_total = architecture_k_of_N(M, K2, pathway_freq, lambda_s_pathway)
                            lambda_mz_total = architecture_k_of_N(M, K2, pathway_freq, lambda_mz_pathway)
                            total_freq = 1-binocdf(K2-1, M, pathway_freq); % probability a pathway is 'affected'
                            %                         h_liab_from_twins_total = twin_concordance_to_heritability( ...
                            %                             lambda_mz_total, lambda_s_total, total_freq, 'ACE')
                            
                            for min_freq = [0.00001 0.0001 0.001 0.01 0.1]
                                for max_freq = [0.05 0.1 0.2 0.4 0.8]
                                    if(max_freq > min_freq)
                                        alpha = max_freq-min_freq; beta = min_freq; % Affine transformation
                                        affine_total_freq = total_freq*alpha+beta;
                                        affine_lambda_s_total = (total_freq^2*alpha^2*lambda_s_total + ...
                                            total_freq*2*alpha*beta+beta^2) / (alpha*total_freq+beta)^2
                                        affine_lambda_mz_total = (total_freq^2*alpha^2*lambda_mz_total + ...
                                            total_freq*2*alpha*beta+beta^2) / (alpha*total_freq+beta)^2
                                        if(affine_lambda_mz_total < affine_lambda_s_total^2)
                                            affine_h_liab_from_twins_total = twin_concordance_to_heritability( ...
                                                affine_lambda_mz_total, affine_lambda_s_total, affine_total_freq, 'ACE')
                                            affine_h_liab_from_mz_twins_total = ...
                                                familial_risk_to_heritability(affine_lambda_mz_total, ...
                                                'liability', affine_total_freq, 1)
                                            affine_h_liab_from_dz_twins_total = ...
                                                familial_risk_to_heritability(affine_lambda_s_total, ...
                                                'liability', affine_total_freq, 0.5)
                                            
                                            if(affine_h_liab_from_twins_total < 1)
                                                p_pathway_tab = compute_k_of_N_joint_tab(2*L, K1, freq);
                                                p_total_tab = compute_k_of_N_joint_tab(M, K2, pathway_freq);
                                                p_x_z_are_one = p_pathway_tab(1,2)*p_total_tab(2,1) / sum(p_total_tab(:,1)) + ...
                                                    p_pathway_tab(2,2)*p_total_tab(2,2) / sum(p_total_tab(:,2));
                                                h_add = 2*N*(p_x_z_are_one - freq*total_freq)^2 / ...
                                                    (freq*(1-freq)*total_freq*(1-total_freq)); % take square correlation
                                                h_liab = heritability_scale_change(h_add, 'liability', total_freq)
                                                percent_over_estimate = ...
                                                    100*(affine_h_liab_from_twins_total - h_liab) / h_liab;
                                            end
                                        end
                                    end
                                end % loop on high-rish penetrance
                            end % loop on baseline penetrance
                        end % loop on number of pathways turned 'on'
                    end % loop on number turned 'on' in a pathway
                end % loop on #genes in pathway
            end % loop on #pathways
        end % loop on frequency
end % switch arch type


