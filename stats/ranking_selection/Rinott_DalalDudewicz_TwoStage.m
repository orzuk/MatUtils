% Analysis of constants in Rinot and Dalal procedures
AssignGeneralConstants;
%two_stage_figs_dir = 'C:\Users\oorzu\Google Drive\Selection procedures\RankingAndSelectionAsymptotics\OperationsResearch\figs';
two_stage_figs_dir = fullfile('C:\Users\', user_str, 'Google Drive\Selection procedures\RankingAndSelectionAsymptotics\RankingSelectionSubmission\EJS\Revision\figs');
line_width = 1.5; % set line width
run_flag=0;

k_vec = unique(round(logspace(log10(1.5), 7, 500))); % 500:500:100000; % list of k values - number of different populations is k+1
pcs_vec = [0.5 0.95]; % [0.5 0.9 0.95 0.99] % [0.1 0.5 0.7 0.9 0.99]; % probability of getting maximum correctly (PCS)
num_k = length(k_vec); % maximum k
nu_vec = [1 2 3 5 10]; % 5 10]; %  100]; % deg. freedom for t-distribution (works for 2 !!!), nu = N0-1

legend_vec = num2str_cell(num2cell(nu_vec));
for i_nu=1:length(nu_vec)
    legend_vec{i_nu} = ['\nu=' legend_vec{i_nu}];
end

p=0.95; nu=1; k=1200;
q_p = (-1/log(p))^(1/nu);
h_1_dalal_approx = (gamma((nu+1)/2) / (nu^(1-nu/2)*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p % Compute approximations
h_1_rinott_approx = 2^(1/nu)*h_1_dalal_approx % Compute approximations
h_1_dalal = fzero(@(x) two_stage_integral_dalal(x, k, nu)-p, [0.3 2] .* h_1_rinott_approx)
h_1_rinott1 = fzero(@(x) two_stage_integral_rinott(x, nu, inf)-(1-p^(1/k)), [0.5 3] .* h_1_dalal_approx)

epsilon=0.01;
x_vec = (tinv(epsilon, nu)-h_1_rinott1):0.1:(-tinv(epsilon, nu));
y_vec = tcdf(x_vec + h_1_rinott1, nu) .* tpdf(x_vec, nu);

if(run_flag)
    [h_k_rinott, h_k_dalal, h_k_rinott_approx, h_k_dalal_approx] = compute_h_k_matrices(pcs_vec, nu_vec, k_vec);    
    % Save run:
    save(fullfile(github_dir, 'MatUtils\stats\ranking_selection\data', 'h1_h2_numerics'), ...
        'h_k_dalal', 'h_k_rinott', 'h_k_dalal_approx', 'h_k_rinott_approx', 'k_vec', 'nu_vec', 'pcs_vec', 'num_k');
else
    load(fullfile(github_dir, 'MatUtils\stats\ranking_selection\data', 'h1_h2_numerics'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: Asymptotics of h_k^j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; rel_error_vec = zeros(num_k, length(pcs_vec)*2); % Figure 1 in paper
[ha, pos] = tight_subplot(3,length(pcs_vec),[.07 .104],[.11 .02],[.1 -.00]); % flip sides !!!
for i_pcs = 1:length(pcs_vec)
    for log_flag=1 % (a.),(b.): show ratio of approximations
        for proc_flag = 0:1 % dalal or rinott
            i_col = 2*i_pcs + proc_flag - 1;
            y_lim = [0 0];
            axes(ha(2*proc_flag+1+(i_pcs-1)));
            for i_nu = 1:length(nu_vec)
                if(proc_flag == 0) % dalal
                    proc_str = 'h_k^{(1)}'; proc_tilde_str = '\tilde{h}_k^{(1)}';
                    rel_error_vec(:,i_col) = reshape((h_k_dalal_approx(i_pcs, i_nu,:)-h_k_dalal(i_pcs, i_nu,:)) ./ h_k_dalal(i_pcs, i_nu,:), num_k, 1);
                else
                    proc_str = 'h_k^{(2)}'; proc_tilde_str = '\tilde{h}_k^{(2)}';
                    rel_error_vec(:,i_col)  = reshape((h_k_rinott_approx(i_pcs, i_nu,:)-h_k_rinott(i_pcs, i_nu,:)) ./ h_k_rinott(i_pcs, i_nu,:), num_k, 1);
                end
                rel_error_table(:,i_col) = interp1(k_vec, rel_error_vec(:,i_col), 10.^(1:7)); 
                y_lim(1) = min(y_lim(1), min(rel_error_vec(:,i_col) ));
                y_lim(2) = max(y_lim(2), max(rel_error_vec(:,i_col) ));
                
                if(log_flag)
                    semilogx(k_vec, rel_error_vec(:,i_col), [color_vec(i_nu)], 'LineWidth', line_width); hold on;
                else
                    plot(k_vec, rel_error_vec(:,i_col), [color_vec(i_nu)], 'LineWidth', line_width); hold on;
                end
            end % loop on i_nu
            xlabel('$k$', 'interpreter', 'latex', 'fontsize', 12); ylabel(['$(' proc_tilde_str '-' proc_str ') / ' proc_str '$'], 'interpreter', 'latex', 'fontsize', 12);
            table_head{i_col} = ['$(' proc_tilde_str '-' proc_str ') / ' proc_str ', p=' num2str(pcs_vec(i_pcs)) '$'];
            if(proc_flag == 0)
                title(['$p=' num2str(pcs_vec(i_pcs)) '$'], 'interpreter', 'latex');
            end
            if((proc_flag == 0) && (i_pcs==1))
%                legend(legend_vec, 'location', 'north', 'fontsize', 7); legend('boxoff');
                legend(legend_vec, 'fontsize', 7, 'position', [0.27 0.82 0.15 0.15]); legend('boxoff');
            end
            y_marg = max(abs(y_lim) * 0.1);
            xlim([0.99 max(k_vec)*1.01]); ylim([y_lim(1) - y_marg, y_lim(2) + y_marg]);
            y_lim = get(gca, 'ylim'); y_lim = min(y_lim, 2); 
            if(i_pcs == length(pcs_vec))
                y_lim = max(min(y_lim, 0.1), -0.1);
            end
            ylim(gca, y_lim);            
            text( 0.9, 0.9, ['(' 'a'+proc_flag+3*(i_pcs-1) '.)'], 'units', 'normalized');
            set(ha(2*proc_flag+1+(i_pcs-1)),'XTick', 10.^(0:7));
            a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
            add_faint_grid(0.7); 
        end % loop on procedure
    end % loop on log
    
    for log_flag=0 % (c.) show ratio
        axes(ha(5+(i_pcs-1)));
        for i_nu = 1:length(nu_vec)
            if(log_flag==0)
                semilogx(k_vec, reshape(((h_k_rinott(i_pcs, i_nu,:) ./ h_k_dalal(i_pcs, i_nu,:))).^2, num_k, 1), ...
                    [color_vec(i_nu) ], 'LineWidth', line_width); hold on;
                ylabel('$(h_k^{(2)} / h_k^{(1)})^2$', 'interpreter', 'latex', 'fontsize', 12);
            else
                semilogx(k_vec, reshape(1./log2((h_k_rinott(i_pcs, i_nu,:) ./ h_k_dalal(i_pcs, i_nu,:)).^2), num_k, 1), ...
                    [color_vec(i_nu) '*'], 'LineWidth', line_width); hold on; %   [color_vec(i_nu)], 'linestyle', symbol_vec{ceil(i_nu/6)} ); hold on;
                ylabel('$1 / 2 \log_2(\frac{h_k^{(2)}}{h_k^{(1)}})$', 'interpreter', 'latex', 'fontsize', 12);
            end
        end % loop on nu
        xlabel('$k$', 'interpreter', 'latex',  'fontsize', 12);        
        if(log_flag==0)
            ylim([1 4]); % [1 2]
        else
            ylim([0.99 max(nu_vec)+1]);
        end
        xlim([0.99 k_vec(end)*1.01]);
        y_lim = get(gca, 'ylim');
        text( 0.9, 0.9, ['(' 'c'+3*(i_pcs-1) '.)'], 'units', 'normalized');
        set(ha(5+(i_pcs-1)),'XTick', 10.^(0:7));
        a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
%         if(i_pcs == length(pcs_vec))
               add_faint_grid(0.7); 
%         end
    end % loop on log flag
end
my_saveas(gcf, fullfile(two_stage_figs_dir, 'h1_and_h2_sqr'), {'epsc', 'jpg'});

relative_error_table_latex = [ [{'$k$ \textbackslash $p$'} table_head]'  [num2cell(10.^(1:7)') num2str_cell(num2cell(rel_error_table), '%2.3f')]' ]';
relative_error_table_latex = latex(relative_error_table_latex, 2)

relative_error_matrix_save = [ [{'Relative error for h_k for two procedures, as function of k and p'} cell(1, 4)]' ...
     [{'p'} num2cell(pcs_vec([1 1 2 2]))]' [{'Proc.', '1', '2', '1', '2'} ]'  [{'k'} cell(1,4)]' [num2cell(k_vec') num2str_cell(num2cell(rel_error_vec), '%2.3f')]' ]';
savecellfile(relative_error_matrix_save, fullfile(github_dir, 'MatUtils\stats\ranking_selection\data', 'h1_h2_numerics.txt'), [], 1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: Optimal choice of nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New: find best h_k
find_best_h=1;  % Figure 3 in paper
if(find_best_h)
    ttt=cputime;
    k_vec = unique(round(logspace(log10(1.5), 6, 20))); % [2 10 100 1000 10000 100000]; % unique(round(logspace(log10(1.5), 7, 500))); % 500:500:100000; % list of k values - number of different populations is k+1
    pcs_vec = [0.5 0.9 0.95 0.99] % [0.1 0.5 0.7 0.9 0.99]; % probability of getting maximum correctly (PCS)
    num_k = length(k_vec); % maximum k
    
    % here optimize h EXACTLY !!
    run_best_h=0;    
    pcs = [0.5 0.9 0.95 0.99]; % probability of getting maximum correctly (PCS)
    num_k = length(k_vec); % maximum k
    % New: show approximation for fixed k as nu goes to inft
    k = [100000]; % unique(round(logspace(log10(1.5), 7, 500))); % 500:500:100000; % list of k values - number of different populations is k+1
    nu_vec = unique(round(logspace(log10(1.5), 7, 500))); % 5 10]; %  100]; % deg. freedom for t-distribution (works for 2 !!!), nu = N0-1
    for i_run = 0% :1
        if(i_run==0)
            Delta=1;
            sigma2_vec = 1;
        else
            Delta = 1/1000;
            sigma2_vec = 1./(1:max(k_vec));
        end
        data_file_name = fullfile(github_dir, 'MatUtils\stats\ranking_selection\data', ['h1_h2_best_h_new' num2str(i_run)]);
        if(run_best_h)
            [best_nu_mat, best_h_mat, best_N_mat, ...
                best_nu_mat_approx, best_h_mat_approx, best_N_mat_approx, best_N_mat_using_nu_approx, ...
                best_nu_mat_approx_semi, best_h_mat_approx_semi, best_N_mat_approx_semi] = ...
                deal(zeros(length(k_vec), length(pcs_vec)));
            for i_k=1:length(k_vec)
                optimize_k = i_k
                if(i_run==1)
                    sigma2_vec = 1./(1:k_vec(i_k));
                end
                for i_pcs=1:length(pcs_vec)
                    p=pcs_vec(i_pcs); k=k_vec(i_k);
                    [best_nu_mat(i_k, i_pcs), best_N_mat(i_k, i_pcs)] = ...
                        fmincon(@(x) two_stage_sample_size(p, x, k, 'dalal', 'numeric', Delta, sigma2_vec), 2*log(k_vec(i_k)), -1, -1);
                    best_h_mat(i_k, i_pcs) = two_stage_compute_h_k(p, best_nu_mat(i_k, i_pcs) , k, 'dalal', 'numeric');
                    %                best_N_mat(i_k, i_pcs) = sum(max(best_nu_mat(i_k, i_pcs)+2, best_h_mat(i_k, i_pcs).^2 .* (sigma_vec .^ 2 ./ Delta.^2)));
                    [best_nu_mat_approx_semi(i_k, i_pcs), best_N_mat_approx_semi(i_k, i_pcs)]  = ...
                        fminsearch(@(x) two_stage_compute_h_k(p, x, k, 'dalal', 'asymptotic'), 2*log(k_vec(i_k)));
                    best_h_mat_approx_semi(i_k, i_pcs) = two_stage_compute_h_k(p, best_nu_mat_approx_semi(i_k, i_pcs) , k, 'dalal', 'asymptotic');
                end
            end
            [h_k_rinott3, h_k_dalal3, h_k_rinott_approx3, h_k_dalal_approx3] = ...
                compute_h_k_matrices(pcs, nu_vec(1:end), k);
            %        save('h1_h2_best_h_new', 'h_k_dalal', 'h_k_rinott', 'h_k_dalal_approx', 'h_k_rinott_approx', ...
            save(data_file_name, 'h_k_dalal3', 'h_k_rinott3', 'h_k_dalal_approx3', 'h_k_rinott_approx3', ...
                'best_nu_mat', 'best_h_mat', 'best_N_mat', ...
                'best_nu_mat_approx', 'best_h_mat_approx', 'best_N_mat_approx', 'best_N_mat_using_nu_approx', ...
                'best_nu_mat_approx_semi', 'best_h_mat_approx_semi', 'best_N_mat_approx_semi', ...
                'k_vec', 'nu_vec', 'pcs_vec', 'sigma2_vec', 'Delta');
        else
            load(data_file_name);
            for i_k=1:length(k_vec) % here can run easy computations
                optimize_k = i_k
                for i_pcs=1:length(pcs_vec)
                    p=pcs_vec(i_pcs); k=k_vec(i_k);
                    best_nu_mat_approx(i_k, i_pcs) = 2*log(k_vec(i_k)); % fminsearch(@(x) two_stage_compute_h_k(p, x, k, 'dalal', 'asymptotic'), 1);
                    best_h_mat_approx(i_k, i_pcs) = sqrt(2*exp(1)*log(k_vec(i_k))); % two_stage_compute_h_k(p, best_nu_mat_approx(i_k, i_pcs) , k, 'dalal', 'asymptotic');
                    best_N_mat_approx(i_k, i_pcs) = sum(max(best_nu_mat_approx(i_k, i_pcs)+2, 2*exp(1)*log(k_vec(i_k)) .* sigma2_vec ./ Delta^2)) .* k_vec(i_k) / length(sigma2_vec); % best_h_mat_approx(i_k, i_pcs).^2 .* k;
                    best_N_mat_using_nu_approx(i_k, i_pcs) = two_stage_sample_size(p, best_nu_mat_approx(i_k, i_pcs), k, 'dalal', 'numeric', Delta, sigma2_vec);
                end
            end
            save(data_file_name, 'h_k_dalal3', 'h_k_rinott3', 'h_k_dalal_approx3', 'h_k_rinott_approx3', ...
                'best_nu_mat', 'best_h_mat', 'best_N_mat', ...
                'best_nu_mat_approx', 'best_h_mat_approx', 'best_N_mat_approx', 'best_N_mat_using_nu_approx', ...
                'best_nu_mat_approx_semi', 'best_h_mat_approx_semi', 'best_N_mat_approx_semi', ...
                'k_vec', 'nu_vec', 'pcs_vec', 'sigma2_vec', 'Delta');
        end
        if(i_run==0)
            figure; % Plot figure: nu. Figure 3 in paper
            [ha, pos] = tight_subplot(1, 3, [.07 .0954],[.11 .02],[.1 -.00]); % flip sides !!! (3,2)
%            [ha, pos] = tight_subplot(3, 1, [.07 .0954],[.11 .02],[.1 -.00]); % flip sides !!! (3,2)
        end
        
        legend_vec = num2str_cell(num2cell(pcs_vec));
        for i_pcs=1:length(pcs_vec)
            legend_vec{i_pcs} = ['$p=' legend_vec{i_pcs} '$'];
        end        
        axes(ha(1+i_run)); %     subplot(1,3,3); % Plot figure: sample size
        for i_pcs = 1:length(pcs)
            loglog(nu_vec, h_k_dalal3(i_pcs,:), [color_vec(i_pcs)], 'linewidth', line_width); hold on;
        end
        loglog(nu_vec, h_k_dalal_approx3(1,:), 'color', orange, 'linestyle', '--', 'linewidth', line_width);  % Approx
        xlabel('$\nu$', 'interpreter', 'latex');  ylabel('$h_k^{(1)}(\nu)$', 'interpreter', 'latex');
        %    legend({'$h_k^1(\nu)$', '$\tilde{h}_k^1(\nu)$'}, 'location', 'northwest', 'interpreter', 'latex'); legend('boxoff');
        for i_pcs = 1:length(pcs)
            I_min = find(h_k_dalal3(i_pcs,:).^2 < nu_vec+2, 1);
            loglog(nu_vec(I_min), h_k_dalal3(i_pcs,I_min), [color_vec(i_pcs) '*'], 'linewidth', line_width);  % Approx
            loglog(nu_vec(I_min), h_k_dalal3(i_pcs,I_min), [color_vec(i_pcs) 'd'], 'linewidth', line_width);
        end
        [~, I_min_approx] = min(h_k_dalal_approx3(1,:));
        loglog(nu_vec(I_min_approx), h_k_dalal_approx3(1,I_min_approx), '*', 'color', orange, 'linewidth', line_width);  % Approx
        loglog(nu_vec(I_min_approx), h_k_dalal_approx3(1,I_min_approx), 'd', 'color', orange, 'linewidth', line_width);  % Approx
        %    I_min2 = find(h_k_dalal3 < nu_vec+2, 1); loglog(nu_vec(I_min), h_k_dalal3(I_min), 'r*', 'linewidth', 2);  % Approx
        xlim([1, max(nu_vec)*1.01]); ylim([3 200]);
        text( 0.84, 0.975, '(a.)', 'units', 'normalized');
        set(ha(1+i_run),'XTick', 10.^(0:7));
        a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
        add_faint_grid(0.7);

        axes(ha(2+i_run)); %  3+  subplot(1,3,1);
        for i_pcs = 1:length(pcs_vec)
            semilogx(k_vec, best_nu_mat(:,i_pcs), [color_vec(i_pcs)], 'linewidth', line_width); hold on;
            % new: fit linear regression:
            beta_mat_for_nu(:,i_pcs) = polyfit(log(k_vec(2:end)'), best_nu_mat(2:end,i_pcs), 1)
        end
        %    title('Optimal $\nu$', 'interpreter', 'latex');
        xlabel('$k$', 'interpreter', 'latex'); ylabel('$\nu_k^{(1*)}$', 'interpreter', 'latex');
        semilogx(k_vec, best_nu_mat_approx(:,i_pcs), 'color', orange, 'linestyle', '--', 'linewidth', line_width); hold on; % optimal nu %  2*log(k_vec)
        should_be_beta_nu = polyfit(log(k_vec'), 2*log(k_vec'), 1)
        xlim([1, max(k_vec)*1.01]); ylim([0 57]); 
        text( 0.84, 0.975, '(b.)', 'units', 'normalized');
        set(ha(2+i_run),'XTick', 10.^(0:7)); %3+
        a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
        add_faint_grid(0.7);

        %    semilogx(k_vec, sqrt(2*exp(1)*log(k_vec)), ['y--'], 'linewidth', 2); hold on; % optimal h_k^1
        
        axes(ha(3+i_run)); % 5+ subplot(1,3,2); % Plot figure: sample size
        for i_pcs = 1:length(pcs_vec)
            semilogx(k_vec, best_N_mat(:,i_pcs)./vec2column(k_vec), [color_vec(i_pcs)], 'linewidth', line_width); hold on;
            beta_mat_for_N(:,i_pcs) = polyfit(log(k_vec(2:end)'), best_N_mat(2:end,i_pcs)./k_vec(2:end)', 1)
        end
        %    title('Optimal $N$', 'interpreter', 'latex');
        xlabel('$k$', 'interpreter', 'latex'); ylabel('$\mu_k^{(1*)}$', 'interpreter', 'latex');
        %    legend(legend_vec, 'location', 'northeast', 'interpreter', 'latex'); legend('boxoff');
        semilogx(k_vec, best_N_mat_approx(:, 1)./vec2column(k_vec), 'color', orange, 'linestyle', '--', 'linewidth', line_width); hold on; % optimal N % 2*exp(1).*log(k_vec)
        semilogx(k_vec, best_N_mat_using_nu_approx(:, 1)./vec2column(k_vec), 'color', orange, 'linestyle', '-', 'linewidth', line_width); hold on; % optimal nu
        should_be_beta_N = polyfit(log(k_vec'), 2*exp(1).*log(k_vec'), 1);
        legend([legend_vec '$\sim$'], 'location', 'northwest', 'interpreter', 'latex'); legend('boxoff');
        xlim([1, max(k_vec)*1.01]);
        text( 0.84, 0.975, '(c.)', 'units', 'normalized');
        set(ha(3+i_run),'XTick', 10.^(0:7)); %5+
        a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
        add_faint_grid(0.7);
        %    semilogx(k_vec, sqrt(2*exp(1)*log(k_vec)), ['y--'], 'linewidth', 2); hold on; % optimal h_k^1
        
    end % loop on i_run
    my_saveas(gcf, fullfile(two_stage_figs_dir, 'h1_and_h2_best'), {'epsc', 'jpg'});
    total_run_time = cputime-ttt
end % find best h


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End figures for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Other figures:
additional_figures=0;
if(additional_figures)
    for(log_flag = 0:1)
        for i_pcs = 1:length(pcs_vec)
            figure;
            for i_nu = 1:length(nu_vec)
                subplot(3, 2, i_nu);
                if(log_flag)
                    loglog(k_vec, h_k_dalal(i_pcs, i_nu,:) , '*'); hold on;
                    loglog(k_vec, h_k_rinott(i_pcs, i_nu,:) , 'r*');
                    loglog(k_vec, h_k_dalal_approx(i_pcs, i_nu,:) , 'b');
                    loglog(k_vec, h_k_rinott_approx(i_pcs, i_nu,:) , 'r');
                else
                    plot(k_vec, h_k_dalal(i_pcs, i_nu,:) , '*'); hold on;
                    plot(k_vec, h_k_rinott(i_pcs, i_nu,:) , 'r*');
                    plot(k_vec, h_k_dalal_approx(i_pcs, i_nu,:) , 'b');
                    plot(k_vec, h_k_rinott_approx(i_pcs, i_nu,:) , 'r');
                end
                xlabel('k'); ylabel('$h_k^i$', 'interpreter', 'latex');
                title(['$\nu=' num2str(nu_vec(i_nu)) '$'], 'interpreter', 'latex');
                
                if(i_nu == length(nu_vec)) % last plot
                    legend({'$h_k^1$', '$h_k^2$', '$h_k^1$ (approx.)', '$h_k^2$ (approx.)'}, 'location', 'southeast', 'interpreter', 'latex'); legend('boxoff');
                end
            end % loop on nu
            %    supertitle('asdfsaf');
            supertitle(['$h_k^1, h_k^2$ for $PCS=' num2str(pcs_vec(i_pcs)) '$'], 'interpreter', 'latex');
            %%        my_saveas(gcf, fullfile(two_stage_figs_dir, ['h1_and_h2_pcs_' strrep(num2str(pcs_vec(i_pcs)), '.', '_')]), {'epsc', 'jpg'});
        end % loop on pcs
    end
    
    
    simulate_two=0;% Now implement procedures and call them:
    if(simulate_two)
        N0=100; PCS = 0.5; Delta = 0.1; n_pop = 50; mu_vec = zeros(1, n_pop); mu_vec(end) =  Delta; sigma_vec = ones(1, n_pop); iters = 100;
        [max_I, PCS_D, h_D] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'dalal');
        [max_I_R, PCS_R, h_R] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'rinott');
        
        PCS_D
        PCS_R
    end
    
    
    
    % Check monotonicity
    nu1=2; nu2=15; h=0.5; k1=1; k2=1; k=14;
    ret = quadgk(@(t) tcdf(t+h, nu1).^k1 .* tcdf(t, nu1).^(k2-1) .* tpdf(t, nu1) - tcdf(t+h, nu2).^k1 .* tcdf(t, nu2).^(k2-1) .* tpdf(t, nu2)  , -inf, inf, 'AbsTol', 10^(-15), 'RelTol', 10^(-10))
    t_vec = -10:0.01:10;
    figure; plot(t_vec, tcdf(t_vec, nu1).^k, 'b', 'linewidth', 2); hold on;
    plot(t_vec, tcdf(t_vec, nu2).^k, 'r', 'linewidth', 2);
    xlim([-0.5 0.5]);
    nu1=4; nu2=4.3; h=0.5; k1=1; k2=1; k=1;
    t_vec2 = 0:0.01:10; figure; y_vec2 = tcdf(h+t_vec2, nu2).^k+tcdf(h-t_vec2, nu2).^k - ...
        tcdf(h+t_vec2, nu1).^k-tcdf(h-t_vec2, nu1).^k;
    plot(t_vec2, y_vec2, 'r', 'linewidth', 2);
    title(['min=' num2str(min(y_vec2))]);
    
    % Show two pdf's - do they really cross at zero?
    nu1=1; nu2=1.2; h=0.5; k1=3; k2=4; k=14;
    h_vec = -10:0.05:10; z_vec1 = zeros(size(h_vec)); z_vec2=z_vec1; z_vec3=z_vec1;
    for i_h = 1:length(h_vec)
        if(mod(i_h, 100) == 0)
            run_i = i_h
        end
        z_vec1(i_h) = two_stage_integral_dalal(h_vec(i_h), [k1 k2], [nu1 nu1]);
        z_vec2(i_h) = two_stage_integral_dalal(h_vec(i_h), [k1 k2], [nu2 nu2]);
        z_vec3(i_h) = two_stage_integral_dalal(h_vec(i_h), [k1 k2], [nu2 nu1]);
    end
    
    %z_vec1 = normcdf(h_vec, 1, 3);
    %z_vec2 = normcdf(h_vec, 1, 2);
    figure; plot(h_vec, z_vec1, 'b', 'linewidth', 2); hold on;
    plot(h_vec, z_vec2, 'r', 'linewidth', 2);
    plot(h_vec, z_vec3, 'g', 'linewidth', 2);
    xlabel('h'); ylabel('$G_{\nu}^k(h)$', 'interpreter', 'latex');
    legend({['\nu_1=' num2str(nu1)], ['\nu_2=' num2str(nu2)], 'mixed'}); legend('boxoff');
    xlim([-0.5 0.5]);
    
    % Look at difference of t random variables:
    nu1=2; nu2=5;  k=1; h=5.5;
    w_vec = tcdf(h+t_vec2, nu2).^k + tcdf(h-t_vec2, nu2).^k - (tcdf(h+t_vec2, nu1).^k + tcdf(h-t_vec2, nu1).^k);
    figure; plot(t_vec2, w_vec);
    figure; plot(t_vec2, tcdf(h+t_vec2, nu2).^k + tcdf(h-t_vec2, nu2).^k, 'b'); hold on;
    plot(t_vec2, tcdf(h+t_vec2, nu1).^k + tcdf(h-t_vec2, nu1).^k, 'r');
    xlabel('t'); ylabel('diff diff');
    title(['min=' num2str(min(w_vec))]);
    figure; plot(t_vec, tcdf(h+t_vec, nu2).^k,  'b'); hold on;
    plot(t_vec, tcdf(h+t_vec, nu1).^k, 'r');
    
    
    % Look at G_nu^k(h+t)+G_nu^k(h-t) as function of t for h>0, t>0
    nu1=1; nu2=5;  k=4; h=4.5;
    u_vec = tcdf(h+t_vec2, nu2).^k + tcdf(h-t_vec2, nu1).^k;
    figure; plot(t_vec2, u_vec);
    xlabel('t'); ylabel('sum');
    title(['min=' num2str(min(u_vec))]);
    
    figure; plot(t_vec, tcdf(t_vec-h, nu1), 'b'); hold on;
    plot(t_vec, tcdf(t_vec-h, nu2), 'b--');
    plot(t_vec, k .* tcdf(t_vec, nu1).^(k-1) .* tpdf(t_vec, nu1), 'r');
    
end % if additional figures

% % %     [h_k_rinott, h_k_dalal, h_k_rinott_approx, h_k_dalal_approx] = compute_h_k_matrices(pcs_vec, nu_vec, k_vec);
% % %     save('h1_h2_best_h', 'h_k_dalal', 'h_k_rinott', 'h_k_dalal_approx', 'h_k_rinott_approx', 'k_vec', 'nu_vec', 'pcs_vec');
% % % %load('h1_h2_best_h');
% % %
% % %
% % %     figure; % Plot figure for
% % %     for i_k = 1:length(k_vec)
% % %         subplot(3,2,i_k);
% % %         for i_pcs = 1:length(pcs_vec)
% % %             semilogy(nu_vec, h_k_dalal(i_pcs,:, i_k), [color_vec(i_pcs) '--'], 'linewidth', 2); hold on;
% % %         end
% % %         for i_pcs = 1:length(pcs_vec)
% % %             semilogy(nu_vec, h_k_dalal_approx(i_pcs,:, i_k), color_vec(i_pcs), 'linewidth', 2); hold on;
% % %         end
% % %         title(['k=' num2str(k_vec(i_k))]);
% % %         xlabel('$\nu$', 'interpreter', 'latex'); ylabel('$h_k^1(\nu)$', 'interpreter', 'latex');
% % %         if(i_k == 1)
% % %             legend(legend_vec, 'location', 'northeast', 'interpreter', 'latex'); legend('boxoff');
% % %         end
% % %     end

