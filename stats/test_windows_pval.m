
% Test the p-values for significance, using both windows methods, 
% and see which one is better
ttt = cputime;

N = 1000; % Size of interval 
n = 140;   % maximal number of points
iters = 5; % number of iterations
min_dist = 2; % minimal distance to compare 
FDR = -1;


order_pvals = zeros(1, iters); 
binom_pvals = zeros(1, iters); 
order_pvals_ten = zeros(1, iters); 
binom_pvals_ten = zeros(1, iters); 
order_pvals_five = zeros(1, iters); 
binom_pvals_five = zeros(1, iters); 

all_order_pvals = zeros(iters, n);
all_binom_pvals = zeros(iters, n);


for i=1:iters
    
    progress = i
    
    for j=n-2:n  % Run over the number of points
        points =  sort(floor(rand(1,j)*N)+1);
        
        temp = set_windows_pval(floor((points+1)/2), N/2, min_dist, FDR); 
        all_order_pvals(i,j) = temp(1,3);
    
% %         temp = simple_binom_windows_pval(points, N, min_dist, FDR); 
% %         all_binom_pvals(i,j) = temp(1,3);
% %                         
        jjj = j
    end
             
end

% Now sort them 
all_order_pvals = sort(all_order_pvals);
all_binom_pvals = sort(all_binom_pvals);

save order_pvals_upto_50_points_file 'all_order_pvals'
save binom_pvals_upto_50_points_file 'all_binom_pvals'
    
figure; hold on; title('order p-vals'); xlabel('points number'); ylabel('p-value'); imagesc(all_order_pvals); colorbar;
% % figure; hold on; title('binom p-vals'); xlabel('points number'); ylabel('p-value'); imagesc(all_binom_pvals); colorbar;

% % % % %     %Draw the points uniformly at random 
% % % % %     points = sort(floor(rand(1,n)*N)+1);
% % % % %     temp = set_windows_pval(points, N, min_dist, FDR); 
% % % % %     order_pvals(i) = temp(1,3);
% % % % %     
% % % % %     temp = simple_binom_windows_pval(points, N, min_dist, FDR); 
% % % % %     binom_pvals(i) = temp(1,3);
% % % % %         
% % % % %     points = sort(floor(rand(1,10)*N)+1);
% % % % %     temp = set_windows_pval(points, N, min_dist, FDR); 
% % % % %     order_pvals_ten(i) = temp(1,3);
% % % % %     
% % % % %     temp = simple_binom_windows_pval(points, N, min_dist, FDR); 
% % % % %     binom_pvals_ten(i) = temp(1,3);
% % % % %     
% % % % %     points = sort(floor(rand(1,5)*N)+1);
% % % % %     temp = set_windows_pval(points, N, min_dist, FDR); 
% % % % %     order_pvals_five(i) = temp(1,3);
% % % % %     
% % % % %     temp = simple_binom_windows_pval(points, N, min_dist, FDR); 
% % % % %     binom_pvals_five(i) = temp(1,3);
% % % % %     
% % % % % 
% % % % %     
% % % % %     
% % % % %     if(mod(i, 20) == 0)
% % % % %         progress = i    
% % % % %     end    
% % % % % end

% % % % % %%order_pvals = binom_pvals;
% % % % % figure; hold on; plot([0:1/(iters-1):1], sort(order_pvals), 'g+');  plot([0:1/(iters-1):1], sort(binom_pvals), 'r*'); 
% % % % % plot([0:1/(iters-1):1], sort(order_pvals_ten), 'mx');  plot([0:1/(iters-1):1], sort(binom_pvals_ten), 'bo');
% % % % % plot([0:1/(iters-1):1], sort(order_pvals_five), 'cs');  plot([0:1/(iters-1):1], sort(binom_pvals_five), 'kd');
% % % % % 
% % % % % title('comparison of two methods p-vals');
% % % % % xlabel('p'); ylabel('pvals'); legend('order-statistics', 'simple-binom', 'order-statistics ten', 'simple-binom ten', 'order-statistics five', 'simple-binom five');
% % % % % 
% % % % % 
% % % % % figure; hold on; plot(log([0:1/(iters-1):1]), log(sort(order_pvals)), 'g+');  plot(log([0:1/(iters-1):1]), log(sort(binom_pvals)), 'r*'); 
% % % % % plot(log([0:1/(iters-1):1]), log(sort(order_pvals_ten)), 'mx');  plot(log([0:1/(iters-1):1]), log(sort(binom_pvals_ten)), 'bo');
% % % % % plot(log([0:1/(iters-1):1]), log(sort(order_pvals_five)), 'cs');  plot(log([0:1/(iters-1):1]), log(sort(binom_pvals_five)), 'kd');
% % % % % 
% % % % % title('LOG LOG comparison of two methods p-vals');
% % % % % xlabel('p'); ylabel('pvals'); legend('order-statistics', 'simple-binom', 'order-statistics ten', 'simple-binom ten', 'order-statistics five', 'simple-binom five');


cputime - ttt