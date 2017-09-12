% Compute the orderering of randomly generated samples 
function F = CheckPermTopFunc(rand_flag, one_side_flag, true_corr_flag, dist_sig, nsamples, ...
    num_vars, alpha, iters, f_res, R, soft_constrain_flag, miu, prior)
 
TRUE = 1; FALSE = 0;            

UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;
TOL = 0.000000000001;
% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0; 
% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0; 


% Here never sample from data. Only from the model ! 
[kept_frac_dist f_all_genes f_all_genes_one_NTOP] = sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, dist_sig, nsamples, num_vars, alpha, iters, R, FALSE);
frac_mean = mean(kept_frac_dist)/(alpha*num_vars);

% Take linear with C_max = 0.4
[inf_limit_frac inf_limit_std x_alpha] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, dist_sig, nsamples, alpha, miu, prior)
%%%5x_alpha = x_alpha * 0.1632993162

f_star = inf_limit_frac; % just another name for convention

% New ! Here we try to compute the saddle point approximation (first order)
% numerically. This is a generalization of the infinite limit fraction
% The last two variables f_star and x_alpha are for numerical technical
% reasons, given as a starting point for the numerical calculations
[saddle_kept_frac_dist sigma_expansion] = ...
    compute_saddle_point_first_order_estimation(rand_flag, one_side_flag, true_corr_flag, dist_sig, ...
    nsamples, num_vars, alpha, f_star, x_alpha,f_res,  soft_constrain_flag); 

this_integral_should_be_one = f_res * sum(saddle_kept_frac_dist)
sigma_we_got_is = sigma_expansion;
f_vec = [f_res:f_res:1-f_res]; 
sigma_numeric = sqrt( sum((saddle_kept_frac_dist ./ sum(saddle_kept_frac_dist)) .* (f_vec.^2)) - ...
                     (sum((saddle_kept_frac_dist ./ sum(saddle_kept_frac_dist)) .* f_vec))^2  )

% Now print the figure. The simulation together with the delta of ngenes->infinity 
if(one_side_flag == ONE_SIDE)
    side_str = 'one sided, ';
else
    side_str = 'two sides, ';
end
if(true_corr_flag == TRUE_AND_SAMPLED)
    corr_str = 'true and sampled, ';
else
    corr_str = 'two sampled, ';
end
if(rand_flag == GAUSSIAN)
    rand_str = 'gaussian corrs, ';
end
if(rand_flag == UNIFORM)
    rand_str = 'uniform corrs, ';
end
if(rand_flag == LINEAR)
    rand_str = 'linear corrs, ';
end

figure; hold on; title(['Kept frac. dist. sim. ' rand_str corr_str side_str 'Ngenes=' num2str(num_vars) ' alpha=' num2str(alpha) ' Nsamp=' num2str(nsamples) ' sim-iters=' num2str(iters)]);
xlabel('f* Kept Frac.'); ylabel('Prob.');
num_bins = min(200, max(kept_frac_dist)-min(kept_frac_dist)+1); % Important ! Added a bin ! 
[hh hh_bins ] = hist(kept_frac_dist./ (alpha*num_vars), num_bins); hh = hh/sum(hh); % Calc the hist and normalize to give probabilities
bin_size = hh_bins(2)-hh_bins(1);
hh_bins = hh_bins + 0.5*(hh_bins(2)-hh_bins(1)); % Adjust bins due to matlab bug !!!! 
bar(hh_bins, hh); %%hist(kept_frac_dist./ (alpha*num_vars), num_bins); 
plot([inf_limit_frac inf_limit_frac], [0 max(hh)*1.2], 'r'); plot([frac_mean frac_mean], [0 max(hh)*1.2], 'g');  
legend('inf. Ngenes approx.', 'empirical mean');

sprintf('Finished first sampling')


% figure;  title ('Here comes the saddle'); xlabel('f_star'); ylabel('P(f_star)');
% plot(saddle_kept_frac_dist, '+');


% Now try a combined plot : Plot only relevant (non-zero) area
start_index = min(find(saddle_kept_frac_dist > TOL));
end_index = max(find(saddle_kept_frac_dist > TOL));
f_vec = [f_res:f_res:1-f_res];
% generate a gaussian curve
points_vec = [f_vec(start_index):(f_vec(end_index)-f_vec(start_index))/1000:f_vec(end_index)];
gauss_points_vec = (1./(sqrt(2.*pi).* sigma_expansion)) .* exp(-(points_vec-f_star).^2./(2.*sigma_expansion.^2));


figure; hold on; title(['Kept frac. Saddle and dist. sim. ' rand_str corr_str side_str 'Ngenes=' num2str(num_vars) ' alpha=' num2str(alpha) ' Nsamp=' num2str(nsamples) ' sim-iters=' num2str(iters)]);
xlabel('f* Kept Frac.'); ylabel('Prob.');

plot(f_vec(start_index:end_index), saddle_kept_frac_dist(start_index:end_index),  '+');  
plot(points_vec, gauss_points_vec, 'c');
plot( hh_bins, hh ./ bin_size, '*r'); 
legend('saddle', 'gaussian expansion', 'simulated');

% Dummy figure not to distroy the previous one !!!! 
figure; 

F = 1; % Dummy assignment
