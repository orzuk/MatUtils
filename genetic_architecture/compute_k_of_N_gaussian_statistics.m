% Compute the correlation between one gaussian and max of many gaussians.
% There is no analytical solution for the general case - the computations
% are done via simulations
%
% Input:
% mu_vec - Gaussians means
% sigma_vec - Gaussians st.d.
% h_x - heritability of each gaussian
% h_share_env - shared envirounment for different gaussians
% iters - # of iterations for computations
% operation - kind of genetic architecture (default is 'MAX')
% operation_param - parameter used for operation 
% N_vec - vector of number of loci (i.e. number of total liabilities) 
% num_bins - for histogram
% compute_mode - how to compute output (default is 'numeric' - numerical integration)
% isoheritability_flag - if this flag is 'ON', we need to make the total heritabilities equal for different N's (not the single Gaussians)
% h_pop_str - what types of heritability population estimates to output 
% 
% Output:
% mu - qtl output mean
% sigma - qtl output st.d.
% k_of_N_hist - qtl histogram probabilities
% k_of_N_bins_loc - qtl histogram values
% corr_vec - correlation of genotype (each gaussian) to phenotype
% h_all - heritability explained by best additive model
% h_pop - cell array of different heritability population estimators 
% qtl_R - correlation in QTL values (phenotypes) between different family members
% same_inds_prob - probability that two relatives have the same Gaussian as maximal
%
function [mu sigma k_of_N_hist k_of_N_bins_loc corr_vec ...
    h_all h_pop h_x_vec qtl_R same_inds_prob] = ...
    compute_k_of_N_gaussian_statistics(mu_vec, sigma_vec, h_x, h_shared_env, iters, ...
    operation, operation_param, ...
    N_vec, num_bins, compute_mode, isoheritability_flag, h_pop_str)


if(~exist('compute_mode', 'var') || isempty(compute_mode))
    compute_mode = 'numeric'; % 'simulation'; % 'numeric'; % compute via numerical integration
end
if(~exist('num_bins', 'var') || isempty(num_bins))
    num_bins = 500;
end
if(~exist('h_pop_str', 'var') || isempty(h_pop_str)) % set default estimator 
    h_pop_str = {'ACE'};
end

N = length(mu_vec); % number of different Gaussians
if(~exist('h_x', 'var') || isempty(h_x))
    h_x = 1;
end
if(~exist('h_shared_env', 'var') || isempty(h_shared_env)) % Default is NO shared environment (need to check that this works!!!) 
    h_shared_env = 0;
end
if(~exist('operation', 'var') || isempty(operation))
    operation = 'MAX'; % operation = 'GxE'; % 'MAX'; %'PROD'; % 'SUM';
end

h_x_vec = zeros(length(N_vec), 1); h_shared_env_vec = zeros(length(N_vec),1); 
k_R_vec = [1 0.5]; % stop at siblings %  0.25 0.125 0.0625 0.03125]; % kinship coefficients
max_family_degree = length(k_R_vec);

if(isoheritability_flag) % here make QTL heritability equal
    h_z = h_x + h_shared_env; % we require TOTAL QTL narrow-sense (epidemiological) heritability to be set
    if(h_z > 1)
        error('Wrong Input Parameters! - total variance cannot exceed one!!!!')
        return;
    end
    for i=1:length(N_vec) % perform binary search seperately for each N
        run_N = N_vec(i)
        % % %         [h_pop qtl_R same_inds_prob] = ... % compute just for twins
        % % %             qtl_familial_correlations_internal(N_vec(i), N_vec(i), k_R_vec(1), operation, ...
        % % %             h_x, h_shared_env, compute_mode, iters); % compute familial correlations and epidemiological heritability
        minbnd_time = cputime;
        % Here assume that BOTH h_x and h_shared_environment are on the
        % Trait's scale (not a single liability scale!)
        switch compute_mode 
            case {'semi-analytic', 'numeric'} % narrow search range 
                lower_bound = h_z;
            otherwise
                lower_bound = 0; 
        end

        tmp_r_MZ_one_liability = fminbnd(@(x) ...
            abs(qtl_familial_correlations_internal(N_vec(i), N_vec(i), k_R_vec(1), operation, ...
            operation_param, x, 0*h_shared_env, compute_mode, iters, 'MZ')-h_z), ...
            lower_bound, 1); % we know that h_z < h_x < 1 but h_pop can be MORE than h_x !!!!
        switch compute_mode
            case {'semi-analytic', 'numeric'}
                upper_bound = tmp_r_MZ_one_liability; 
                lower_bound = h_x;
            otherwise
                upper_bound = 1; 
                lower_bound = 0;
        end
        h_x_vec(i) = fminbnd(@(x) ...
            abs(qtl_familial_correlations_internal(N_vec(i), N_vec(i), k_R_vec(1), operation, ...
            operation_param, x, 0*h_shared_env, compute_mode, iters, 'MZ')-(h_x + 0*h_shared_env)), ...
            lower_bound, upper_bound); % we know that h_z < h_x < 1 but h_pop can be MORE than h_x !!!!
        
%        h_x_vec(i) = 2*(tmp_r_MZ_one_liability-tmp_r_DZ_one_liability); 
        h_shared_env_vec(i) = max(0,tmp_r_MZ_one_liability-h_x_vec(i)); % tmp_r_DZ_one_liability; 
        minbnd_time = cputime - minbnd_time
               
        %         xxx = 0.5; ttt = cputime;
        %         qtl_familial_correlations_internal(N_vec(i), N_vec(i), k_R_vec(1:2), operation, ...
        %             xxx, h_shared_env, compute_mode, iters)
        %         cputime - ttt        
    end
else
    h_x_vec(:) = h_x; h_shared_env_vec(:) = h_shared_env;
end

block_size = 50000; % don't allow more simulations at once
num_blocks = iters / min(iters, block_size); % assumed integer

same_inds_prob = zeros(max_family_degree,length(N_vec));
corr_vec = zeros(N,length(N_vec)); % compute correlation of each gaussian to phenotype (QTL)
x_vec = -10:0.001:10;

switch compute_mode % compute heritability via simulation or numeric integration 
    case {'simulation', 'simulations'}        
        sum_data_genetic = zeros(1, N); sum_sqr_data_genetic = zeros(1,N);
        sum_qtl_vec = zeros(1,length(N_vec)); sum_sqr_qtl_vec = zeros(1,length(N_vec));
        sum_data_genetic_times_qtl_vec = zeros(N, length(N_vec));
    case {'semi-analytic', 'numeric'}
        corr_vec = zeros(N,length(N_vec)); % compute correlation of each gaussian to phenotype (QTL)
        Exz = zeros(length(N_vec),1); 
        mu = zeros(length(N_vec),1); sigma = zeros(length(N_vec),1); 
        h_all = zeros(length(N_vec),1); 
        for i=1:length(N_vec)  % Compute integration numerically (no sampling is needed). Works only for a few architectures 
            switch operation
                case {'MIN', 'MAX'}
                    corr_vec(i) = quadl(@(x)max_and_single_gaussian_integrand_internal( ...
                        x, N_vec(i)), -5, 10); % Here goes the actual computation
                    Exz(i) = corr_vec(i);
                    [mu(i) sigma(i)] = maxnormstat(N_vec(i)); % compute moments via numeric integrtion (no sampling).
                    if(strcmp(operation, 'MIN')) % flip order
                        mu(i) = -mu(i);
                    end
                case 'DIFF' % here z = |x2-x1|^p
                    corr_vec(i) = 0; % correlation of x1 and |x2-x1|^p: E x1*|x2-x1|^p
                    Exz(i) = corr_vec(i);
                    [mu(i) sigma(i)] = powernormstat(operation_param(1)); % compmute moments of z
            end % switch operation
            corr_vec(i) = corr_vec(i)*sqrt(h_x_vec(i)) / sigma(i); % normalize
            h_all(i) = N_vec(i)*corr_vec(i)^2; % compute heritability (take into account liability)            
        end % loop on N_vec         
end

debug_figures = 0;

switch compute_mode
    case {'simulation', 'simulations'}
        for b=1:num_blocks
            data_genetic = randn(block_size, N) .* sqrt(h_x) .* ...
                repmat(sigma_vec, block_size, 1) + repmat(mu_vec, block_size, 1); % simulate data. (different Gaussians)
            data = data_genetic + randn(block_size, N) .* ...
                sqrt(1-h_x) .* repmat(sigma_vec, block_size, 1);  % environmental component of each gaussian !
            qtl_vec = qtl_operation_internal(data, operation, operation_param, N_vec);  % compute qtl vec for ALL N's from 1 to N
            
            sum_data_genetic = sum_data_genetic + sum(data_genetic);
            sum_sqr_data_genetic = sum_sqr_data_genetic + sum(data_genetic.^2);
            sum_qtl_vec = sum_qtl_vec + sum(qtl_vec,2)';
            sum_sqr_qtl_vec = sum_sqr_qtl_vec + sum(qtl_vec .^2,2)';
            for i=1:N % compute empirical correlations
                for j=1:length(N_vec)
                    if(i <= N_vec(j))
                        sum_data_genetic_times_qtl_vec(i,j) = sum_data_genetic_times_qtl_vec(i,j) + ...
                            sum(data_genetic(:,i) .* qtl_vec(j,:)');
                        
                    end
                end
            end
            
            if(debug_figures) % Perform a different calculation to get correlation vec. Works currently for max of identical Gaussians.
                x = data(:,1); % take one Gaussian
                for i=1:length(N_vec)
                    if(N_vec(i) > 1)
                        y_vec = (N_vec(i)-1) .* normcdf(x_vec).^(N_vec(i)-2) .* normpdf(x_vec); % Randomize the maximum of N-1 Gaussians
                        y = distrnd(x_vec, y_vec, block_size, 1); % compute y: the maximum of N-1 Gaussians
                        corr_vec2(i) = corr(x, max(x,y)) * sqrt(h_x); % compute correlation
                    else
                        corr_vec2(i) = sqrt(h_x);
                    end
                    h_all3(i) = corr_vec2(i)^2*N_vec(i); % compute heritability
                end
            end % debug figures
            for i=1:length(N_vec)
                [~,~,~,~,STATS] = regress(vec2column(qtl_vec(i,:)), ...
                    [data_genetic(:,1:N_vec(i)) ones(block_size, 1)]); % compute heritability via regression
                h_all2(i) = STATS(1);
            end
        end % loop on blocks        
        mean_data_genetic = sum_data_genetic ./ iters;
        std_data_genetic = sqrt(sum_sqr_data_genetic ./ iters - mean_data_genetic.^2);
        mean_qtl_vec = sum_qtl_vec ./ iters;
        std_qtl_vec = sqrt(sum_sqr_qtl_vec ./ iters - mean_qtl_vec.^2);
        mean_data_genetic_times_qtl_vec = sum_data_genetic_times_qtl_vec ./ iters;
                        
        for i=1:N % compute empirical correlations
            for j=1:length(N_vec)
                if(i <= N_vec(j))
                    corr_vec(i,j) = (mean_data_genetic_times_qtl_vec(i,j) - ...
                        mean_data_genetic(i).*mean_qtl_vec(j)) ./ ...
                        (std_data_genetic(i) .* std_qtl_vec(j));
                    %            corr(data_genetic(:,i), vec2column(qtl_vec(j,:))); % correlate only with genetic component!
                end
            end
        end
end % simulate data



switch compute_mode
    case {'simulation', 'simulations'}
        h_all = sum(corr_vec.^2); % take sum of individual loci heritabilities
    case {'semi-analytic', 'numeric'}
        %        h_all = h_all4;
end


if(debug_figures) % Debug figures: compute heritability in different ways:
    figure; hold on; plot(h_all, h_all2, '.');
    plot(h_all, h_all3, 'g.'); plot(h_all, h_all4, 'm.'); plot(0:1,0:1, 'r');
    title('two narrow-sense heritability computations');
    xlabel('single-locus sum'); ylabel('other computations');
    legend({'regression', 'one locus and sum-of-N-1', 'semi-analytic integration: one locus and sum-of-N-1'});
    
    % Compare mean and st.d. estimates
    figure; subplot(2,1,1); hold on; plot(mean(qtl_vec,2), mu4, '.');
    plot(0:3,0:3, 'r'); title('Compare max-qtls means');
    xlabel('\mu empirical'); ylabel('\mu semi-analytic integration');
    subplot(2,1,2); hold on; plot(std(qtl_vec,[],2), sigma4, '.');
    plot(0:3,0:3, 'r'); title('Compare max-qtls st.d.');
    xlabel('\sigma empirical'); ylabel('\sigma semi-analytic integration');
end



k_of_N_hist = cell(length(N_vec),1);
k_of_N_bins_loc = cell(length(N_vec),1);
switch compute_mode
    case {'simulation', 'simulations'}
        mu = mean(qtl_vec,2); sigma = std(qtl_vec,[],2); % compute first two moments
        normalize_flag = 0;
        if(normalize_flag) % normalize phenotype
            qtl_vec = qtl_vec - mean(qtl_vec); qtl_vec = qtl_vec ./ std(qtl_vec);
        end
        for i=1:length(N_vec)
            [k_of_N_hist{i} k_of_N_bins_loc{i}] = hist(qtl_vec(i,:), num_bins); % compute histogram
        end
    case  {'semi-analytic', 'numeric'}
        for i=1:length(N_vec)
            k_of_N_bins_loc{i} = linspace(-10,10,num_bins);
            switch operation
                case 'MIN'
                    k_of_N_hist{i} = maxnormpdf(-k_of_N_bins_loc{i}, N_vec(i)); % take negative (assume standard Gaussians)
                otherwise % we know only how to write max analytically 
                    k_of_N_hist{i} = maxnormpdf(k_of_N_bins_loc{i}, N_vec(i));
            end
                    
        end
end



if(debug_figures)
    %figure; hist_density(qtl_vec, num_bins, 'r');
    %title(['Dist. of QTL values: ' operation ' of Gaussians']);
    
    figure; hold on; plot(log(N_vec), log(Exz), '.');
    plot(log(N_vec), -log(N_vec), 'r')
    title('heritability explained as function of N');
    xlabel('N (log)'); ylabel('Exz  (log)'); % 'h_{loci}');
    
    %figure; hold on; plot(
    
    
    figure; plot(N_vec ./ log(N_vec), h_all ./ Exz, '.'); title('heritability explained as function of N');
    xlabel('N '); ylabel('h_{loci} / Exz'); % 'h_{loci}');
    
    % figure; hold on; (-log(2.*log(N_vec(3:end))./N_vec(3:end)), log(h_all(3:end)), '.'); title('heritability explained as function of N');
    % plot(0:0.1:5, -(0:0.1:5), 'r');
    % xlabel('2log(N)/N'); ylabel('h_{loci}'); legend('h_all (numeric integration)', 'y=x');
    
    
    
    full_figure(0); loglog(N_vec, h_all, '.');  hold on; title('Heritability explained for max-of-N-Gaussians');
    loglog(N_vec, 2*log(N_vec)./N_vec, 'r');
    xlabel('N (log)'); ylabel('h_{loci} (log)');
    legend('exact (numeric integration)', 'asymptotic: 2log N/N');
    my_saveas(gcf, 'heritability_explained_max_of_Gaussians', format_fig_vec);
end % more temporary figs


[~, qtl_R same_inds_prob] = ... % Here h_pop is based on 2*(r_MZ-r_DZ) !!!
    qtl_familial_correlations_internal(N_vec, N, k_R_vec, operation, operation_param, ...
    h_x_vec, h_shared_env_vec', compute_mode, iters, 'ACE'); % compute familial correlations and epidemiological heritability
h_pop = cell(length(h_pop_str),1); 
for i=1:length(h_pop_str) % compute all heritability estimates
    h_pop{i}  = ... % Here h_pop is based on 2*r_DZ-r_MZ !!!
        qtl_familial_correlations_internal(N_vec, N, k_R_vec, operation, operation_param, ...
        h_x_vec, h_shared_env_vec', compute_mode, iters, h_pop_str{i}); % compute familial correlations and epidemiological heritability
end

mu = vec2row(mu); sigma = vec2row(sigma);
h_all = vec2row(h_all); h_pop = vec2row(h_pop); h_x_vec = vec2row(h_x_vec);






% Internal function: Compute maximum and single gaussian correlations. 
% g_N(x) * x (x PHI(x) / (N-1) - phi(x))
function ret = max_and_single_gaussian_integrand_internal(x, N)

%ret = x^2 .* ( normpdf(x) .* normcdf(x).^N + normcdf(x) .* N .* normcdf(x).^(N-1) .* normpdf(x)
ret = maxnormpdf(x,N-1) .* x .* (x.*normcdf(x)./(N-1) - normpdf(x)); % Simplified


switch N
    case 1
        ret(:) = normpdf(x) .* x.^2;
        %    case 2
        %        ret = normpdf(x) .* normcdf(x) .* x.^2 ; % Simplified
        
end





% % % % Cumulative for Debug
% % % function ret = joint_two_gaussian_maxima_cumulative_internal(z, z_R, N, rho)
% % %
% % %
% % % m = size(z_R, 1); n = size(z_R, 2);
% % %
% % % temp_int = zeros(1,n);
% % % for i=1:n % compute integral one by one
% % %     temp_int(i) = quadl(@(x)  norm_density_times_cumulative_internal(x, z_R(i), rho) , -10, z(i));        % That's the slowest part ...
% % %     temp_int2(i) = mvncdf([z(i) z_R(i)]', [0 0]', [1 rho; rho 1]);
% % % end
% % % ret2 = repmat(temp_int.^(N-1), size(z_R, 1), 1);
% % % z_R_normalized = (z_R-z.*rho)./sqrt(1-rho.^2); % normalize z_R conditioned on z
% % % ret2 = ret2 .* normcdf(z_R_normalized);
% % % ret = maxnormcdf(z,N) .* ret2; % multiplication by z*z_R should keep symmetry



