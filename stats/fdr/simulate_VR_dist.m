% Simulates p-values, apply FDR mand compute V and R distributions.
% Calculate statistics of the number of rejections (R), false
% rejections (V) and their ratio (V/R^+)
% 
% Input: 
% iters - number of simulations
% m - number of hypothesis
% m0 - number of null hypothesis
% mu_gap - difference between null and non-null hypothesis
% dep_flag - what kind of dependence between pvals with employ
% rho - correlation coefficients between (Z-transformed) p-vals
% fdr_str - which fdr procedure to use
% q - parameter for fdr procedure
% two_side_flag - parameter for p-values simulations (do one-sided or two-sided test)
% 
% Output: 
% VR_hist - histogram of V and R (2-d)
% V_moments - mean and std of V
% R_moments - mean and std of R
% FDR_moments - mean and std of FDR=V/R^+
%
function [VR_hist V_moments R_moments FDR_moments] = ...
    simulate_VR_dist(iters,m,m0,mu_gap,dep_flag,rho,fdr_str,q, two_side_flag, varargin)

if(exist('two_side_flag', 'var'))
    pvals_vec = simulate_pvals(iters,m,m0,mu_gap,dep_flag,rho, two_side_flag); % simulate p-vals
else
    pvals_vec = simulate_pvals(iters,m,m0,mu_gap,dep_flag,rho); % simulate p-vals
end
if(length(pvals_vec) <200000)
    figure; hold on; subplot(2,2,1); plot(pvals_vec(:,1), pvals_vec(:,end), '.'); xlabel('1st pval'); ylabel('2nd pval');  % plot first vs. last (can be null vs. non-null)
    subplot(2,2,3); hist(pvals_vec(:,1), 500); title('pval 1');
    subplot(2,2,4); hist(pvals_vec(:,2), 500); title('pval 2');
end
 
[R V] = FDR_mat_main(pvals_vec, q, fdr_str, [], m0); % assume no auxillary u for now


V_moments = [mean(V) var(V)]; R_moments = [mean(R) var(R)]; FDR_moments = [mean(V./max(R,1)), var(V./max(R,1))];  

VR_hist = hist2d ([R V], -0.5:1:max(R)+1, -0.5:1:max(R)+1) ./ iters;
fig_flag = 0;
if(fig_flag)
%    figure; Plot2dHist(VR_hist, 0:max(V)+1, 0:max(R)+1, 'V', 'R', '(V,R) two-dim histogram');     % not so clear plot
    figure;
    if(length(VR_hist) == 1)
        imagesc(VR_hist);
    else
        imagesc(VR_hist, [min(VR_hist(:)) max(max(VR_hist(2:end,:)))]);
    end
    xlabel('V(+1)'); ylabel('R(+1)'); title(['(V,R) two-dim histogram. P(V=0,R=0) is actually ' num2str(VR_hist(1,1))]); colorbar;
end

    
