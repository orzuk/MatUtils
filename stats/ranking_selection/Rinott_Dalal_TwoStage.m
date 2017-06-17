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
%h_1_dalal_almost = fsolve(@(x) tcdf(x, nu)^k-p, h_1_dalal_approx)
h_1_rinott1 = fzero(@(x) two_stage_integral_rinott(x, nu, inf)-(1-p^(1/k)), [0.5 3] .* h_1_dalal_approx)


debug_integral=0;
if(debug_integral)
    h_vec = 0:100:10^5; %(h_1_rinott_approx*2);
    int_rinott_vec = zeros(size(h_vec));
    for i_h=1:length(h_vec)
        run_h = h_vec(i_h)
        int_rinott_vec(i_h) = two_stage_integral_rinott(h_vec(i_h), nu, inf)-(1-p^(1/k));
    end
    figure; loglog(h_vec, abs(int_rinott_vec)); hold on;
    
    figure; plot(h_vec, int_rinott_vec); hold on;
    plot(h_1_rinott_approx, 0, 'r*');
    plot(h_1_dalal_approx, 0, '*b');
    ylim([-0.0003 0.0003]);
    xlim([55000 65000]);
    
    two_stage_integral_rinott(60050, nu, [-inf 0])
    two_stage_integral_rinott(60050, nu, [0 inf])
    two_stage_integral_rinott(61050, nu, [-inf 0])
    two_stage_integral_rinott(61050, nu, [0 inf])
    two_stage_integral_rinott(61250, nu, [-inf 0])
    two_stage_integral_rinott(61250, nu, [0 inf])
    two_stage_integral_rinott(62200, nu, inf)
end % debug integral


epsilon=0.01;
x_vec = (tinv(epsilon, nu)-h_1_rinott1):0.1:(-tinv(epsilon, nu));
y_vec = tcdf(x_vec + h_1_rinott1, nu) .* tpdf(x_vec, nu);
figure; plot(x_vec, y_vec); hold on;

%h_1_rinott22 = fzero(@(x) two_stage_integral_rinott(x, nu, inf)-p, [0.5*h_1_dalal_approx, 2*h_1_dalal_approx])
%h_1_rinott2 = fsolve(@(x) two_stage_integral_rinott(x, nu, 20)-p, h_1_dalal_approx)
%h_1_rinott3 = fsolve(@(x) two_stage_integral_rinott(x, nu, -tinv(10^(-7), nu))-p, h_1_dalal_approx)

h_k_rinott = []; h_k_dalal = []; h_k_rinott_approx = []; h_k_dalal_approx = [];
for i_pcs = 1:length(pcs_vec)
    for i_nu = 1:length(nu_vec)
        p = pcs_vec(i_pcs), nu = nu_vec(i_nu)
        q_p = (-1/log(p))^(1/nu); % gevinv(p, nu, 1, 0); % the p-th quantile of nu-Frechet distribution with c.d.f.: e^{-x^{-nu}}
        [h_k_rinott{i_pcs, i_nu}, h_k_rinott_approx{i_pcs, i_nu}, h_k_dalal{i_pcs, i_nu}, h_k_dalal_approx{i_pcs, i_nu}] = deal(zeros(1, num_k+1));
        
        for i_k=1:length(k_vec)
            k=k_vec(i_k);
            if(nu == Inf)
                h_k_dalal_approx{i_pcs, i_nu}(i_k+1) = sqrt(2*log(k));
                h_k_rinott_approx{i_pcs, i_nu}(i_k+1) = 2*sqrt(log(k));
            else
                h_k_dalal_approx{i_pcs, i_nu}(i_k+1) = (gamma((nu+1)/2) / (nu^(1-nu/2)*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p; % Compute approximations
                h_k_rinott_approx{i_pcs, i_nu}(i_k+1) = (gamma((nu+1)/2) / (nu^(1-nu/2)*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p * 2^(1/nu); % sqrt(2); % Compute approximations
            end
            h_k_dalal{i_pcs, i_nu}(i_k+1) = fzero(@(x) two_stage_integral_dalal(x, k, nu)-p, ...
                [-0.01 3] .* h_k_rinott_approx{i_pcs, i_nu}(i_k+1));
            h_k_rinott{i_pcs, i_nu}(i_k+1) = fzero(@(x) two_stage_integral_rinott(x, nu)-(1-p^(1/k)), ... %   p, ... % , [tinv(epsilon, nu)-x -tinv(epsilon, nu)])-p, ...
                [-0.01 3] .* h_k_dalal_approx{i_pcs, i_nu}(i_k+1));
        end
        h_k_dalal{i_pcs, i_nu} = h_k_dalal{i_pcs, i_nu}(2:end);
        h_k_rinott{i_pcs, i_nu} = h_k_rinott{i_pcs, i_nu}(2:end);
        h_k_dalal_approx{i_pcs, i_nu} = h_k_dalal_approx{i_pcs, i_nu}(2:end);
        h_k_rinott_approx{i_pcs, i_nu} = h_k_rinott_approx{i_pcs, i_nu}(2:end);
        
        %figure; plot(k_vec, h_k_dalal, '*');
        %hold on; plot(k_vec, h_k_rinott, 'r*');
        %legend({'dalal', 'rinott'}); legend('boxoff');
        %xlabel('k'); ylabel('h(k)');
    end % loop on nu
end % loop on pcs

% Save run: 
%save('h1_h2_numerics', 'h_k_dalal', 'h_k_rinott', 'h_k_dalal_approx', 'h_k_rinott_approx', 'k_vec', 'nu_vec', 'pcs_vec'); 
%load('h1_h2_numerics'); 


for(log_flag = 0:1)
    for i_pcs = 1:length(pcs_vec)
        figure;
        for i_nu = 1:length(nu_vec)
            subplot(3, 2, i_nu);
            if(log_flag)
                loglog(k_vec, h_k_dalal{i_pcs, i_nu} , '*'); hold on;
                loglog(k_vec, h_k_rinott{i_pcs, i_nu} , 'r*');
                loglog(k_vec, h_k_dalal_approx{i_pcs, i_nu} , 'b');
                loglog(k_vec, h_k_rinott_approx{i_pcs, i_nu} , 'r');
            else
                plot(k_vec, h_k_dalal{i_pcs, i_nu} , '*'); hold on;
                plot(k_vec, h_k_rinott{i_pcs, i_nu} , 'r*');
                plot(k_vec, h_k_dalal_approx{i_pcs, i_nu} , 'b');
                plot(k_vec, h_k_rinott_approx{i_pcs, i_nu} , 'r');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[ha, pos] = tight_subplot(3,length(pcs_vec),[.07 .104],[.11 .02],[.1 -.00]); % flip sides !!!
for i_pcs = 1:length(pcs_vec)
    %    for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
    %set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
    %         set(ha(1:3), 'fontsize', 8);
    
    for log_flag=1 % (a.),(b.): show ratio of approximations
        for proc_flag = 0:1 % dalal or rinott
            y_lim = [0 0];
            axes(ha(2*proc_flag+1+(i_pcs-1))); %% subplot(1, 3, proc_flag+1); %            figure;
            for i_nu = 1:length(nu_vec)
                if(proc_flag == 0) % dalal
                    proc_str = 'h_k^1'; proc_tilde_str = '\tilde{h}_k^1';
                    rel_error_vec = (h_k_dalal_approx{i_pcs, i_nu}-h_k_dalal{i_pcs, i_nu}) ./ h_k_dalal{i_pcs, i_nu};
                else
                    proc_str = 'h_k^2'; proc_tilde_str = '\tilde{h}_k^2';
                    rel_error_vec = (h_k_rinott_approx{i_pcs, i_nu}-h_k_rinott{i_pcs, i_nu}) ./ h_k_rinott{i_pcs, i_nu};
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
            %ylabel(['$\frac{\bar{' proc_str '}-' proc_str '}{' proc_str '}$'], 'interpreter', 'latex');
            %            title(['Asymptotic accuracy for $' proc_str '$, $PCS=' num2str(pcs_vec(i_pcs)) '$'], 'interpreter', 'latex');
            
            if(proc_flag == 0)
                title(['$p=' num2str(pcs_vec(i_pcs)) '$'], 'interpreter', 'latex');
            end
            if((proc_flag == 0) && (i_pcs==1))
                legend(legend_vec, 'location', 'northeast', 'fontsize', 7); legend('boxoff');
            end
            y_marg = max(abs(y_lim) * 0.1);
            xlim([0.99 max(k_vec)*1.01]); ylim([y_lim(1) - y_marg, y_lim(2) + y_marg]);
            y_lim = get(gca, 'ylim');
            %            text( max(k_vec)*0.2, y_lim(2)*0.9-0.05, ['(' 'a'+proc_flag+3*(i_pcs-1) '.)']); % , 'units', 'normalized');
            text( 0.9, 0.9, ['(' 'a'+proc_flag+3*(i_pcs-1) '.)'], 'units', 'normalized');
            
            set(ha(2*proc_flag+1+(i_pcs-1)),'XTick', 10.^(0:7));
            a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
        end
    end
    
    for log_flag=0 % (c.) show ratio
        axes(ha(5+(i_pcs-1))); %% subplot(1, 3, proc_flag+1); %            figure; subplot(1, 3, 3); % one figure for each PCS
        for i_nu = 1:length(nu_vec)
            if(log_flag==0)
                semilogx(k_vec, (h_k_rinott{i_pcs, i_nu} ./ h_k_dalal{i_pcs, i_nu})  , [color_vec(i_nu) ], 'LineWidth', 2); hold on;
                ylabel('$h_k^2 / h_k^1$', 'interpreter', 'latex', 'fontsize', 12);
                %%%        , 'linestyle', symbol_vec{ceil(i_nu/6)} ); hold on;
            else
                semilogx(k_vec, 1./log2(h_k_rinott{i_pcs, i_nu} ./ h_k_dalal{i_pcs, i_nu})  , [color_vec(i_nu) '*'], 'LineWidth', 2); hold on; %   [color_vec(i_nu)], 'linestyle', symbol_vec{ceil(i_nu/6)} ); hold on;
                ylabel('$1 / \log_2(\frac{h_k^2}{h_k^1})$', 'interpreter', 'latex', 'fontsize', 12);
            end
        end % loop on nu
        xlabel('$k$', 'interpreter', 'latex',  'fontsize', 12);
        
        %        title(['Relative efficiency for $PCS=' num2str(pcs_vec(i_pcs)) '$'], 'interpreter', 'latex');
        %        legend(legend_vec, 'location', 'northwest'); legend('boxoff');
        %    plot(k_vec, repmat(sqrt(2), length(k_vec), 1), 'k--');
        if(log_flag==0)
            ylim([1 2]); % [0 2]
        else
            ylim([0.99 max(nu_vec)+1]);
        end
        xlim([0.99 k_vec(end)*1.01]);
        y_lim = get(gca, 'ylim');
        %        text( max(k_vec)*0.2, y_lim(2)*0.95, ['(' 'c'+3*(i_pcs-1) '.)']);
        text( 0.9, 0.9, ['(' 'c'+3*(i_pcs-1) '.)'], 'units', 'normalized');
        set(ha(5+(i_pcs-1)),'XTick', 10.^(0:7));
        a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8);
        
        
        %    my_saveas(gcf, fullfile(two_stage_figs_dir, ['h1_h2_ratio_pcs_' strrep(num2str(pcs_vec(i_pcs)), '.', '_')]), {'epsc', 'jpg'});
        %        my_saveas(gcf, fullfile(two_stage_figs_dir, ['h1_h2_pcs_' strrep(num2str(pcs_vec(i_pcs)), '.', '_')]), {'epsc', 'jpg'});
    end
end
my_saveas(gcf, fullfile(two_stage_figs_dir, 'h1_and_h2'), {'epsc', 'jpg'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End figure for paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



simulate_two=0;% Now implement procedures and call them:
if(simulate_two)
    N0=100; PCS = 0.5; Delta = 0.1; n_pop = 50; mu_vec = zeros(1, n_pop); mu_vec(end) =  Delta; sigma_vec = ones(1, n_pop); iters = 100;
    [max_I, PCS_D, h_D] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'dalal');
    [max_I_R, PCS_R, h_R] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'rinott');
    
    PCS_D
    PCS_R
end


check_t=0;
if(check_t) % Check sum of T student
    m = 1000000; nu=2;
    T = trnd(nu, m, 2);
    T_sum = sum(T, 2);
    k=50;
    quantile(T_sum ./ sqrt(2), 1/k)
    quantile(T(:,1), 1/k)
    
    x_vec = -100:0.1:100;
    figure; hist_density(T_sum ./ sqrt(2), x_vec); hold on;
    hist_density(T(:,1), x_vec, 'r'); hold on;
    
    figure; plot(sort(T_sum ./ sqrt(2)), (1:m)./m); hold on;
    plot(sort(T_sum ./ 2), (1:m)./m, 'g');
    plot(sort(T(:,1)), (1:m)./m, 'r--');
    xlabel('t'); ylabel('G(t)');
    xlim([quantile(T_sum, 1/200), -quantile(T_sum, 1/200)]);
end






