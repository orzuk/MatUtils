% Analysis of constants in Rinot and Dalal procedures 
two_stage_figs_dir = 'C:\Users\oorzu\Google Drive\Selection procedures\RankingAndSelectionAsymptotics\OperationsResearch\figs'; 
AssignGeneralConstants; 

k_vec = 1:100; % list of k values - number of different populations is k+1 
pcs_vec = [0.1 0.5 0.7 0.9 0.99]; % probability of getting maximum correctly (PCS)
num_k = max(k_vec); % maximum k 
nu_vec = [1 2 10 100]; % deg. freedom for t-distribution (works for 2 !!!), nu = N0-1

legend_vec = num2str_cell(num2cell(nu_vec)); 
for i_nu=1:length(nu_vec)
    legend_vec{i_nu} = ['\nu=' legend_vec{i_nu}]; 
end

h_k_rinott = []; h_k_dalal = []; h_k_rinott_approx = []; h_k_dalal_approx = []; 
for i_pcs = 1:length(pcs_vec)
    for i_nu = 1:length(nu_vec)
        p = pcs_vec(i_pcs); nu = nu_vec(i_nu); 
        q_p = (-1/log(p))^(1/nu); % gevinv(p, nu, 1, 0); % the p-th quantile of nu-Frechet distribution with c.d.f.: e^{-x^{-nu}} 
        [h_k_rinott{i_pcs, i_nu}, h_k_rinott_approx{i_pcs, i_nu}, h_k_dalal{i_pcs, i_nu}, h_k_dalal_approx{i_pcs, i_nu}] = deal(zeros(1, num_k+1));
        
        for k=1:length(k_vec)
            h_k_dalal{i_pcs, i_nu}(k+1) = fsolve(@(x) two_stage_integral_dalal(x, k, nu)-p, h_k_dalal{i_pcs, i_nu}(k));  
            h_k_rinott{i_pcs, i_nu}(k+1) = fsolve(@(x) two_stage_integral_rinott(x, k, nu)-p, h_k_rinott{i_pcs, i_nu}(k));
            h_k_dalal_approx{i_pcs, i_nu}(k+1) = (gamma((nu+1)/2) / (nu*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p; % Compute approximations
            h_k_rinott_approx{i_pcs, i_nu}(k+1) = (gamma((nu+1)/2) / (nu*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p * sqrt(2); % Compute approximations        
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
        
for i_pcs = 1:length(pcs_vec)
    figure; 
    for i_nu = 1:length(nu_vec)
        subplot(2, 2, i_nu);  
        loglog(k_vec, h_k_dalal{i_pcs, i_nu} , '*'); hold on; 
        loglog(k_vec, h_k_rinott{i_pcs, i_nu} , 'r*');
        loglog(k_vec, h_k_dalal_approx{i_pcs, i_nu} , 'b');  
        loglog(k_vec, h_k_rinott_approx{i_pcs, i_nu} , 'r');
        legend({'dalal', 'rinott', 'dalal (appr.)', 'rinott (appr.)'}); legend('boxoff');
        xlabel('k'); ylabel('h(k)');
        title(['$h_k^1, h_k^2$ for $PCS=' num2str(pcs_vec(i_pcs)) ', \nu=' num2str(nu_vec(i_nu)) '$'], 'interpreter', 'latex');
    end % loop on nu
    my_saveas(gcf, fullfile(two_stage_figs_dir, ['h1_and_h2_pcs_' num2str(pcs_vec(i_pcs))]), {'epsc', 'jpg'});
end % loop on pcs
  
for i_pcs = 1:length(pcs_vec)
    figure; % one figure for each PCS
    for i_nu = 1:length(nu_vec)
        plot(k_vec, h_k_rinott{i_pcs, i_nu} ./ h_k_dalal{i_pcs, i_nu} , [color_vec(i_nu) '*']); hold on;
    end % loop on nu    
    xlabel('k'); ylabel('h_{rinnot}(k) / h_{dalal}(k)');
    title(['Relative efficiency for PCS=' num2str(pcs_vec(i_pcs))]);
    legend(legend_vec); legend('boxoff');
    plot(k_vec, repmat(sqrt(2), length(k_vec), 1), 'k--');
    ylim([1 2]); 
    my_saveas(gcf, fullfile(two_stage_figs_dir, ['h1_h2_ratio_pcs_' num2str(pcs_vec(i_pcs))]), {'epsc', 'jpg'});
end % loop on pcs

% Now implement procedures and call them: 

N0=100; PCS = 0.5; Delta = 0.1; n_pop = 50; mu_vec = zeros(1, n_pop); mu_vec(end) =  Delta; sigma_vec = ones(1, n_pop); iters = 100; 
[max_I, PCS_D, h_D] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'dalal');
[max_I_R, PCS_R, h_R] = two_stage_selection_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'rinott');
 
PCS_D
PCS_R

