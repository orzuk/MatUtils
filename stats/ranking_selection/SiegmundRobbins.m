%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New: start Siegmound and Robbins:
AssignGeneralConstants;
two_stage_figs_dir = fullfile('C:\Users\', user_str, 'Google Drive\Selection procedures\RankingAndSelectionAsymptotics\RankingSelectionSubmission\EJS\Revision\figs');
line_width = 1.5; % set default
run_flag=0;
if(run_flag) % run calculations (long time)
    k_vec = unique(round(logspace(log10(1.5), 7, 1000))); % 500:500:100000; % list of k values - number of different populations is k+1
    pcs_vec = [0.5 0.9 0.95 0.99]; num_p = length(pcs_vec);  % [0.5 0.9 0.95 0.99] % [0.1 0.5 0.7 0.9 0.99]; % probability of getting maximum correctly (PCS)
    num_k = length(k_vec); % maximum k
    Delta = 1;
    sigma_vec = 1;
    options = optimset('TolX',10^(-19));
    
    alpha_vec = [0 0.25 0.5 1]; beta_vec = [0.5 0.5 0.5 0.5];
    N_approx_vec = cell(length(alpha_vec), 1); N_star_vec = cell(length(alpha_vec), 1);
    for i_alpha=1:length(alpha_vec)
        s_vec = ceil(beta_vec(i_alpha) .* k_vec .^ alpha_vec(i_alpha)); % select s out of k
        N_star_vec{i_alpha} = ones(num_k, num_p);
        C = alpha_vec(i_alpha); % log(s_vec(end)) / log(k_vec(end)); % limit
        N_approx_vec{i_alpha} = 2 .* sigma_vec(1)^2.*(1+sqrt(C)).^2 .* log(k_vec) ./Delta^2; %  log(log(k_vec)); % get asymptotic result
        
        for i_pcs = 1:length(pcs_vec) % loop on p
            for k=1:length(k_vec)
                if(mod(k, 10) == 0)
                    run_k = k
                end
                if(Gaussian_conv_max(Delta * sqrt(1) / sigma_vec, k_vec(k), s_vec(k)) > pcs_vec(i_pcs))
                    N_star_vec{i_alpha}(k, i_pcs) = 1;
                else
                    N_star_vec{i_alpha}(k, i_pcs) = fzero(@(x) Gaussian_conv_max(Delta * sqrt(x) / sigma_vec, k_vec(k), s_vec(k)) - pcs_vec(i_pcs), ...
                        [1, 10*N_approx_vec{i_alpha}(k)], options);
                end
            end
        end % loop on pcs
    end    
    save(fullfile(github_dir, 'MatUtils\stats\ranking_selection\data', 'N_star_numerics'), 'N_star_vec', 'N_approx_vec', 'k_vec', 'alpha_vec', 'beta_vec', 'pcs_vec'); % save data 
else % load results
    load(fullfile(github_dir, 'MatUtils\stats\ranking_selection\data', 'N_star_numerics'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a.) Plot results: figure showing asymptotics for different values of k^alpha
figure;
legend_vec = num2str_cell(num2cell(1./alpha_vec));
for i_alpha=1:length(alpha_vec)
    legend_vec{i_alpha} = ['$s_k=\frac{1}{2} k^{\frac{1}{' legend_vec{i_alpha} '}}$'];
end
legend_vec{1} = '$s_k=1$'; legend_vec{end} = '$s_k=\frac{1}{2} k$';
for log_flag=0 % :1
    subplot(1, 2, 1);
    for i_alpha=1:length(alpha_vec)
        if(log_flag)
            loglog(k_vec, N_approx_vec{i_alpha}, color_vec(i_alpha), 'linewidth', line_width); hold on; % plot approximation according to generalized reult
        else
            semilogx(k_vec, N_approx_vec{i_alpha}, color_vec(i_alpha), 'linewidth', line_width); hold on; % plot approximation according to generalized reult
        end
    end
    for i_alpha=1:length(alpha_vec)
        for i_pcs = 3 % 1:length(pcs_vec) % 1: % show only 0.95
            if(log_flag)
                loglog(k_vec, N_star_vec{i_alpha}(:, i_pcs), [color_vec(i_alpha) '--'], 'linewidth', line_width); hold on; % plot actual PCS
            else
                semilogx(k_vec, N_star_vec{i_alpha}(:, i_pcs), [color_vec(i_alpha) '--'], 'linewidth', line_width); hold on; % plot actual PCS
            end
        end
    end % loop on alpha
    xlabel('$k$', 'interpreter', 'latex'); ylabel('$N_k^*(p), \: \tilde{N}_k^*$ ', 'interpreter', 'latex'); %
    xlim([0.99 max(k_vec)*1.01]);
    legend(legend_vec, 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 10); legend('boxoff');
    set(gca,'XTick', 10.^(0:7));
    text( 0.89, 0.97, '(a.)', 'units', 'normalized');
    a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
    add_faint_grid(0.7);
end % loop on log-flag
% (b.) Figure showing the ratio for different values of p:
subplot(1, 2, 2);
legend_vec = num2str_cell(num2cell(pcs_vec));
for i_pcs=1:length(pcs_vec)
    legend_vec{i_pcs} = ['$p=' legend_vec{i_pcs} '$'];
end
log_flag=1;
for i_alpha = 3 % take sqrt
    for i_pcs = 1:length(pcs_vec) % 1: % show only 0.95
        relative_error_mat(:,i_pcs) = (N_approx_vec{i_alpha}'-N_star_vec{i_alpha}(:, i_pcs)) ./ N_star_vec{i_alpha}(:, i_pcs);
        semilogx(k_vec, relative_error_mat(:,i_pcs), color_vec(i_pcs), 'linewidth', line_width); hold on; % plot approximation according to generalized result
        relative_error_table(:,i_pcs) = interp1(k_vec, relative_error_mat(:,i_pcs)', 10.^(1:7));
    end
end
xlabel('$k$', 'interpreter', 'latex'); ylabel('$\frac{\tilde{N}_k^*-N_k^*(p)}{N_k^*(p)}$', 'interpreter', 'latex'); %
xlim([0.99 max(k_vec)*1.01]); ylim([-1 3.5]);
legend(legend_vec, 'interpreter', 'latex', 'location', 'north', 'fontsize', 10); legend('boxoff');
set(gca,'XTick', 10.^(0:7));
text( 0.89, 0.97, '(b.)', 'units', 'normalized');
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
add_faint_grid(0.7);

my_saveas(gcf, fullfile(two_stage_figs_dir, 'N_star_k_alpha_diff_PCS_new'), {'epsc', 'jpg'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End figure for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Table for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relative_error_table_latex = [ [{'$k$ \textbackslash $p$'} num2cell(pcs_vec)]'  [num2cell(10.^(1:7)') num2str_cell(num2cell(relative_error_table), '%2.3f')]' ]';
relative_error_table_latex = latex(relative_error_table_latex, 2)

relative_error_mat_save = [ [{'Relative error as function of k and p'} cell(1, 4)]'  [{'k \textbackslash p'} num2cell(pcs_vec)]'  [num2cell(k_vec') num2str_cell(num2cell(relative_error_mat), '%2.3f')]' ]';
savecellfile(relative_error_mat_save, fullfile(github_dir, 'MatUtils\stats\ranking_selection\data', 'N_star_numerics.txt'), [], 1); 



