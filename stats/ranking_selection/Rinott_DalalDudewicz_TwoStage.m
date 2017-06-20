% Analysis of constants in Rinot and Dalal procedures
two_stage_figs_dir = 'C:\Users\oorzu\Google Drive\Selection procedures\RankingAndSelectionAsymptotics\OperationsResearch\figs';
AssignGeneralConstants;

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
figure; plot(x_vec, y_vec); hold on;

[h_k_rinott, h_k_dalal, h_k_rinott_approx, h_k_dalal_approx] = compute_h_k_matrices(pcs_vec, nu_vec, k_vec); 

% Save run: 
save('h1_h2_numerics', 'h_k_dalal', 'h_k_rinott', 'h_k_dalal_approx', 'h_k_rinott_approx', 'k_vec', 'nu_vec', 'pcs_vec'); 
%load('h1_h2_numerics'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[ha, pos] = tight_subplot(3,length(pcs_vec),[.07 .104],[.11 .02],[.1 -.00]); % flip sides !!!
for i_pcs = 1:length(pcs_vec)
    for log_flag=1 % (a.),(b.): show ratio of approximations
        for proc_flag = 0:1 % dalal or rinott
            y_lim = [0 0];
            axes(ha(2*proc_flag+1+(i_pcs-1))); 
            for i_nu = 1:length(nu_vec)
                if(proc_flag == 0) % dalal
                    proc_str = 'h_k^1'; proc_tilde_str = '\tilde{h}_k^1';
                    rel_error_vec = reshape((h_k_dalal_approx(i_pcs, i_nu,:)-h_k_dalal(i_pcs, i_nu,:)) ./ h_k_dalal(i_pcs, i_nu,:), num_k, 1);
                else
                    proc_str = 'h_k^2'; proc_tilde_str = '\tilde{h}_k^2';
                    rel_error_vec = reshape((h_k_rinott_approx(i_pcs, i_nu,:)-h_k_rinott(i_pcs, i_nu,:)) ./ h_k_rinott(i_pcs, i_nu,:), num_k, 1);
                end
                y_lim(1) = min(y_lim(1), min(rel_error_vec));
                y_lim(2) = max(y_lim(2), max(rel_error_vec));
                
                if(log_flag)
                    semilogx(k_vec, rel_error_vec , [color_vec(i_nu)], 'LineWidth', 2); hold on;
                else
                    plot(k_vec, rel_error_vec , [color_vec(i_nu)], 'LineWidth', 2); hold on;
                end
            end % loop on i_nu
            xlabel('$k$', 'interpreter', 'latex', 'fontsize', 12); ylabel(['$(' proc_tilde_str '-' proc_str ') / ' proc_str '$'], 'interpreter', 'latex', 'fontsize', 12);
            if(proc_flag == 0)
                title(['$p=' num2str(pcs_vec(i_pcs)) '$'], 'interpreter', 'latex');
            end
            if((proc_flag == 0) && (i_pcs==1))
                legend(legend_vec, 'location', 'northeast', 'fontsize', 7); legend('boxoff');
            end
            y_marg = max(abs(y_lim) * 0.1);
            xlim([0.99 max(k_vec)*1.01]); ylim([y_lim(1) - y_marg, y_lim(2) + y_marg]);
            y_lim = get(gca, 'ylim');
            text( 0.9, 0.9, ['(' 'a'+proc_flag+3*(i_pcs-1) '.)'], 'units', 'normalized');
            set(ha(2*proc_flag+1+(i_pcs-1)),'XTick', 10.^(0:7));
            a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
        end
    end
    
    for log_flag=0 % (c.) show ratio
        axes(ha(5+(i_pcs-1))); 
        for i_nu = 1:length(nu_vec)
            if(log_flag==0)
                semilogx(k_vec, reshape(((h_k_rinott(i_pcs, i_nu,:) ./ h_k_dalal(i_pcs, i_nu,:))).^2, num_k, 1), ...
                    [color_vec(i_nu) ], 'LineWidth', 2); hold on;
                ylabel('$(h_k^2 / h_k^1)^2$', 'interpreter', 'latex', 'fontsize', 12);
            else
                semilogx(k_vec, reshape(1./log2((h_k_rinott(i_pcs, i_nu,:) ./ h_k_dalal(i_pcs, i_nu,:)).^2), num_k, 1), ...
                    [color_vec(i_nu) '*'], 'LineWidth', 2); hold on; %   [color_vec(i_nu)], 'linestyle', symbol_vec{ceil(i_nu/6)} ); hold on;
                ylabel('$1 / 2 \log_2(\frac{h_k^2}{h_k^1})$', 'interpreter', 'latex', 'fontsize', 12);
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
    end
end
my_saveas(gcf, fullfile(two_stage_figs_dir, 'h1_and_h2_sqr'), {'epsc', 'jpg'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End figure for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Other figures: 
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


% New: find best h_k
find_best_h=1;
if(find_best_h)
    k_vec = [2 10 100 1000 10000 100000]; % unique(round(logspace(log10(1.5), 7, 500))); % 500:500:100000; % list of k values - number of different populations is k+1
    pcs_vec = [0.5 0.9 0.95 0.99] % [0.1 0.5 0.7 0.9 0.99]; % probability of getting maximum correctly (PCS)
    num_k = length(k_vec); % maximum k
    nu_vec = [1:100]; % 5 10]; %  100]; % deg. freedom for t-distribution (works for 2 !!!), nu = N0-1
    
    [h_k_rinott, h_k_dalal, h_k_rinott_approx, h_k_dalal_approx] = compute_h_k_matrices(pcs_vec, nu_vec, k_vec); 
    
    
    figure; % Plot figure
    legend_vec = num2str_cell(num2cell(pcs_vec));
    for i_pcs=1:length(pcs_vec)
        legend_vec{i_pcs} = ['$p=' legend_vec{i_pcs} '$'];
    end
    for i_k = 1:length(k_vec)
        subplot(3,2,i_k);
        for i_pcs = 1:length(pcs_vec)
            semilogy(nu_vec, h_k_dalal(i_pcs,:, i_k), [color_vec(i_pcs) '--'], 'linewidth', 2); hold on;
        end
        for i_pcs = 1:length(pcs_vec)
            semilogy(nu_vec, h_k_dalal_approx(i_pcs,:, i_k), color_vec(i_pcs), 'linewidth', 2); hold on;
        end
        title(['k=' num2str(k_vec(i_k))]);
        xlabel('$\nu$', 'interpreter', 'latex'); ylabel('$h_k^1(\nu)$', 'interpreter', 'latex');
        if(i_k == 1)
            legend(legend_vec, 'location', 'northeast', 'interpreter', 'latex'); legend('boxoff');
        end
    end
    
    % here optimize h
    best_nu_mat = zeros(length(k_vec), length(pcs_vec));
    for i_k=1:length(k_vec)
        optimize_k = i_k
        for i_pcs=1:length(pcs_vec)
            p=pcs_vec(i_pcs); k=k_vec(i_k);
            best_nu_mat(i_k, i_pcs) = fminsearch(@(x) (gamma((x+1)/2) / (x^(1-x/2)*gamma(x/2)*sqrt(pi))) ^ (1/x) * k^(1/x) * (-1/log(p))^(1/x), 1);
        end
    end
    figure; % Plot figure
    for i_pcs = 1:length(pcs_vec)
        semilogx(k_vec, best_nu_mat(:,i_pcs), [color_vec(i_pcs)], 'linewidth', 2); hold on;
    end
    title('Optimal $\nu$', 'interpreter', 'latex');
    xlabel('$k$', 'interpreter', 'latex'); ylabel('$\nu$', 'interpreter', 'latex');
    legend(legend_vec, 'location', 'northeast', 'interpreter', 'latex'); legend('boxoff');
    semilogx(k_vec, 2*log(k_vec), ['y--'], 'linewidth', 2); hold on; % optimal nu
%    semilogx(k_vec, sqrt(2*exp(1)*log(k_vec)), ['y--'], 'linewidth', 2); hold on; % optimal h_k^1
end


simulate_two=0;% Now implement procedures and call them:
if(simulate_two)
    N0=100; PCS = 0.5; Delta = 0.1; n_pop = 50; mu_vec = zeros(1, n_pop); mu_vec(end) =  Delta; sigma_vec = ones(1, n_pop); iters = 100;
    [max_I, PCS_D, h_D] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'dalal');
    [max_I_R, PCS_R, h_R] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'rinott');
    
    PCS_D
    PCS_R
end



