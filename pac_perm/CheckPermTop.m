% Test the order of randomly generated samples. 
% We use here both simulations and analytic expansions. 

path(path, 'E:\Research\PACPerm\numeric\nmm\linalg');
path(path, 'E:\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric');
path(path, 'C:\Weizmann\Research\PACPerm\numeric\nmm\linalg');

% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;
rand_flag = GAUSSIAN; 
TOL = 0.00000000001;

% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0; 
one_side_flag = TWO_SIDES;   % Currently we work only with gaussian two-sided

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0; 
true_corr_flag  = TRUE_AND_SAMPLED;

ttt = cputime;

top_frac_vec = [0.001,0.005,0.025,0.1];  % Vector of alpha fraction of total genes
sigma = 0.5; % the noise per one sampling

max_nsamples = 400; % Maximal number of samplings. We assume that the st.d. is proportional to 1/sqrt(samples)

res = 20; % the resulution of N that we use

nsamples_vec = res:res:max_nsamples;

times_to_sample = 50; % How many times to sample for each value of N


%%%%% Maple Integral 
%%%% int( (1/sqrt(2*Pi)) *  exp(-c*c/2) * int( (1/(sigma*sqrt(2*Pi))) * exp(-x*x/(2*sigma*sigma)),x=x_alpha-c..infinity), c=c_alpha..infinity);

C_alpha = norminv(1-top_frac_vec);  % Get the C_alpha vector

% Save a value for each alpha and each N 
maple_int_res = zeros(length(top_frac_vec),length(nsamples_vec));

sigma_vec = sigma ./ sqrt(nsamples_vec); % st.d. is proportional to 1/sqrt(N) 

for i = 1:length(nsamples_vec)
    for j=1:length(top_frac_vec)
        %sigma_vec(i)
        %C_alpha(j)
        %            F = inline('(1/sqrt(2.*pi)) .* exp(-c.*c./2) .* (1 - normcdf( (C_alpha(j).*(1+sigma)-c)./sigma_vec(i) ))');
        maple_int_res(j,i) = (1.0/top_frac_vec(j))*quad('JointDensFrac', C_alpha(j), 99999, [], [], sigma_vec(i), C_alpha(j), one_side_flag);
    end
end


% One means replace the delta function by an exponential soft constrain.
% Zero means use delta integral representation
soft_constrain_flag = 0; 
alpha = 0.3; num_vars = 1000; nsamples=1; iters = 1000; f_res = 0.01;
% Generate a sampling of the distribution for a given set of parameters 

miu = []; prior = []; % Dummy for now
CheckPermTopFunc(rand_flag, one_side_flag, true_corr_flag, 1, nsamples, num_vars, alpha, iters, f_res, soft_constrain_flag, miu, prior);


% % % [kept_frac_dist f_all_genes]= sample_kept_fraction_distribution(rand_flag, one_side_flag, true_corr_flag, 1, nsamples, num_vars, alpha, iters, R, FALSE);
% % % frac_mean = mean(kept_frac_dist)/(alpha*num_vars);
% % % 
% % % % Take linear with C_max = 0.4
% % % [inf_limit_frac inf_limit_std x_alpha] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, 1, nsamples, alpha)
% % % %%%5x_alpha = x_alpha * 0.1632993162
% % % 
% % % f_star = inf_limit_frac; % just another name for convention
% % % 
% % % % New ! Here we try to compute the saddle point approximation (first order)
% % % % numerically. This is a generalization of the infinite limit fraction
% % % % The last two variables f_star and x_alpha are for numerical technical
% % % % reasons, given as a starting point for the numerical calculations
% % % saddle_kept_frac_dist = compute_saddle_point_first_order_estimation(rand_flag, one_side_flag, true_corr_flag, 1, nsamples, num_vars, alpha, f_star, x_alpha,f_res, soft_constrain_flag); 
% % % 
% % % this_integral_should_be_one = f_res * sum(saddle_kept_frac_dist)
% % % 
% % % % Now print the figure. The simulation together with the delta of ngenes->infinity 
% % % if(one_side_flag)
% % %     side_str = 'one sided, ';
% % % else
% % %     side_str = 'two sides, ';
% % % end
% % % if(true_corr_flag)
% % %     corr_str = 'true and sampled, ';
% % % else
% % %     corr_str = 'two sampled, ';
% % % end
% % % if(rand_flag == GAUSSIAN)
% % %     rand_str = 'gaussian corrs, ';
% % % end
% % % if(rand_flag == UNIFORM)
% % %     rand_str = 'uniform corrs, ';
% % % end
% % % if(rand_flag == LINEAR)
% % %     rand_str = 'linear corrs, ';
% % % end
% % % 
% % % figure; hold on; title(['Kept frac. dist. sim. ' rand_str corr_str side_str 'Ngenes=' num2str(num_vars) ' alpha=' num2str(alpha) ' Nsamp=' num2str(nsamples) ' sim-iters=' num2str(iters)]);
% % % xlabel('f* Kept Frac.'); ylabel('Prob.');
% % % num_bins = min(200, max(kept_frac_dist)-min(kept_frac_dist)+1); % Important ! Added a bin ! 
% % % [hh hh_bins ] = hist(kept_frac_dist./ (alpha*num_vars), num_bins); hh = hh/sum(hh); % Calc the hist and normalize to give probabilities
% % % bin_size = hh_bins(2)-hh_bins(1);
% % % hh_bins = hh_bins + 0.5*(hh_bins(2)-hh_bins(1)); % Adjust bins due to matlab bug !!!! 
% % % bar(hh_bins, hh); %%hist(kept_frac_dist./ (alpha*num_vars), num_bins); 
% % % plot([inf_limit_frac inf_limit_frac], [0 max(hh)*1.2], 'r'); plot([frac_mean frac_mean], [0 max(hh)*1.2], 'g');  
% % % legend('inf. Ngenes approx.', 'ampirical mean');
% % % 
% % % sprintf('Finished first sampling')
% % % 
% % % 
% % % figure;  title ('Here comes the saddle'); xlabel('f_star'); ylabel('P(f_star)');
% % % plot(saddle_kept_frac_dist, '+');
% % % 
% % % 
% % % % Now try a combined plot : Plot only relevant (non-zero) area
% % % start_index = min(find(saddle_kept_frac_dist > TOL));
% % % end_index = max(find(saddle_kept_frac_dist > TOL));
% % % figure; hold on; title(['Kept frac. Saddle and dist. sim. ' rand_str corr_str side_str 'Ngenes=' num2str(num_vars) ' alpha=' num2str(alpha) ' Nsamp=' num2str(nsamples) ' sim-iters=' num2str(iters)]);
% % % xlabel('f* Kept Frac.'); ylabel('Prob.');
% % % f_vec = [f_res:f_res:1-f_res];
% % % plot(f_vec(start_index:end_index), saddle_kept_frac_dist(start_index:end_index),  '+');  plot( hh_bins, hh ./ bin_size, '*r'); legend('saddle','simulated');

total_time = cputime - ttt

return;

% This is the sampling part !!! Should be transferd into a function !!!!!! 
for num_vars = [1000, 10000] %, 10000] 
    %%%%top_vars_vec = [1:10:50:100:150:300:500:1000]; % How many genes are in the top set 
    %%%top_vars_vec = [10,50,200,1000]; % How many genes are in the top set 
    
    top_vars_vec = num_vars .* top_frac_vec;
    
    mean_flipped = zeros(1, length(nsamples_vec)); std_flipped = zeros(1, length(nsamples_vec));
    mean_top_frac = zeros(length(top_vars_vec), length(nsamples_vec)); std_top_frac = zeros(length(top_vars_vec), length(nsamples_vec));
    samptime = 0; sorttime = 0; fliptime = 0; 

    
    for i=1:times_to_sample
        i        
        
        % Sample the true correlations
        corr_vec = sample_true_correlations(rand_flag, num_vars, 0, 1); % Take zero mean and one variance
        
        sss = cputime; 
        
        % Generate the sampling correlation vector
        sample_corr_vec = sample_noisy_correlations(corr_vec, nsamples_vec, 1*GAUSSIAN); % This function should return the sample correlation
        
        samptime = samptime + cputime-sss;
        % Now we need to check various things (e.g. if the perm is perfect, how
        % many errors etc.)
        
        sss = cputime; 
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 Start
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Flip
        % 
        %     % count fraction of flipped pairs
        %     curr_count_flipped = count_flipped_pairs(sample_corr_vec);        
        %     mean_flipped = mean_flipped  + curr_count_flipped';
        %     std_flipped = std_flipped  + (curr_count_flipped').^2;
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 End
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Flip        
        
        fliptime = fliptime + cputime-sss;
        
        sss = cputime; 

        % count fraction of top genes kept
        [ val_corr ind_corr] = sort(sample_corr_vec, 2);
        sorttime = sorttime + cputime-sss;
        
        % Find out how many of the top genes stayed top 
        j=1;
        for top_vars = top_vars_vec
            curr_top_frac = (sum(ind_corr(:, num_vars-top_vars+1:end) > num_vars-top_vars, 2))';
            mean_top_frac(j,:) = mean_top_frac(j,:) +  curr_top_frac;
            std_top_frac(j,:) = std_top_frac(j,:) +  curr_top_frac.^2;
            j=j+1;
        end                
        
    end
    
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 Start
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Flip
    % mean_flipped = 2*mean_flipped ./ ( num_vars * (num_vars-1) .* times_to_sample);
    % std_flipped = sqrt( 4*std_flipped ./ ( num_vars^2 * (num_vars-1)^2 .* times_to_sample) - mean_flipped.^2 );
    % figure; hold on; errorbar(nsamples_vec, mean_flipped, std_flipped); title('Total fraction of flipped corrs'); xlabel('N'); ylabel('Fraction');
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 End
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Flip
    
    j=1;
    for top_vars = top_vars_vec
        mean_top_frac(j,:) = mean_top_frac(j,:) ./ ( top_vars .* times_to_sample);
        std_top_frac(j,:) = sqrt( std_top_frac(j,:) ./ ( top_vars^2 .* times_to_sample) - mean_top_frac(j,:).^2 );
        j=j+1;
    end
    
    
    figure; hold on; 
    
    legend(num2str(top_vars_vec(1)),  num2str(top_vars_vec(2)),  num2str(top_vars_vec(3)), num2str(top_vars_vec(4)));  
    
    % %     errorbar(nsamples_vec, mean_top_frac(1,:), std_top_frac(1,:)); 
    % %     errorbar(nsamples_vec, mean_top_frac(2,:), std_top_frac(2,:), 'r'); 
    % %     errorbar(nsamples_vec, mean_top_frac(3,:), std_top_frac(3,:), 'm'); 
    % %     errorbar(nsamples_vec, mean_top_frac(4,:), std_top_frac(4,:), 'c'); 
    plot(nsamples_vec, mean_top_frac(1,:), ':'); 
    plot(nsamples_vec, mean_top_frac(2,:), 'r:'); 
    plot(nsamples_vec, mean_top_frac(3,:), 'm:'); 
    plot(nsamples_vec, mean_top_frac(4,:), 'c:'); 
    
    
    
    % Now plot the same for the maples
    plot(nsamples_vec, maple_int_res(1,:));
    plot(nsamples_vec, maple_int_res(2,:), 'r');
    plot(nsamples_vec, maple_int_res(3,:), 'm');
    plot(nsamples_vec, maple_int_res(4,:), 'c');
    
    title(['Total fraction of top genes out of ' num2str(num_vars) ' total genes. Dashed is simulation, Solid, solution for infinite number of genes']); xlabel('N'); ylabel('Fraction');
    legend(num2str(top_vars_vec(1)),  num2str(top_vars_vec(2)),  num2str(top_vars_vec(3)), num2str(top_vars_vec(4)));  
    
    samptime_is = samptime
    fliptime_is = fliptime
    sorttime_is = sorttime
    cputime - ttt
    
end

