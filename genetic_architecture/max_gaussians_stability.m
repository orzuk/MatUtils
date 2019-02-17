% Check best embryo (different project) 
r2_vec = 0:0.01:1; n=10; C_alpha = norminv(1-1/n);  TOL = 0.00000000001; M = 1000000;
for j=1:length(r2_vec)
    r2=r2_vec(j); sigma = sqrt(2*(1-r2)/r2);
    r2_is = r2
    f_star_approx(j) = n*quadl('JointDensFrac', C_alpha, 100*(1+sqrt(sigma)),TOL, [], sigma, C_alpha, 1);
    X = mvnrnd(zeros(n,1), eye(n)*r2/2, M); E = mvnrnd(zeros(n,1), eye(n)*(1-r2), M);
    [~, I] = max(X, [], 2); [~, J] = max(X+E, [], 2); f_star_sim(j) = mean(I == J);
end
f_star_approx(1) = 1/n; % fix NaNs
figure; hold on; 
plot(r2_vec, f_star_sim, '-', 'linewidth', 2); 
plot(r2_vec, f_star_approx, 'r-', 'linewidth', 2);
legend({'sim.', 'approx.'}, 'location', 'southeast'); legend('boxoff'); xlabel('r^2'); ylabel('Pr(max-chosen)'); 
add_faint_grid(0.6); 
my_saveas(gcf, ['prob_max_embryo_n=' num2str(n)], {'epsc', 'pdf', 'jpg'}); 




