% Internal function: compute familial correlations - this also gives an
% epidemiological estimation for the heritability
%
% Input:
% N_vec - vector of different N values
% N - maximum N value
% k_R_vec - IBD coefficient
% operation - architecture type (currently only MAX of Gaussians is supported)
% operation_param - paremetrizing architecture 
% h_x_vec - heritability of each liability
% h_shared_env - shared envirounment for sibs/MZ on each liability ! 
% compute_mode - how to compute: 'simulations' or 'semi-analyitc'/'numeric'
% iters - number of iterations when computing using simulations
% epi_str - what epidemiological model to use to estimate heritability
%           (default: 'MZ'. Can also be 'ACE', 'ADE', ..)
% 
% Output:
% h_pop - estimation of heritability (from r_MZ^2)
% qtl_R - table of familial trait correlations
% same_inds_prob - probability that the same liability was chosen as maximum
%
function [h_pop qtl_R same_inds_prob] = ...
    qtl_familial_correlations_internal(N_vec, N, k_R_vec, operation, operation_param, ...
    h_x_vec, h_shared_env, compute_mode, iters, epi_str)

if(~exist('epi_str', 'var') || isempty(epi_str))
    epi_str = 'MZ'; % how to estimate epidemiological heritability 
end
max_family_degree = length(k_R_vec);
qtl_R = zeros(max_family_degree,length(N_vec));  % qtl_R(i,j) is the correlation coefficient of relatives of degree i when N=N_vec(j)
%qtl_R2 = qtl_R; % This one is computed via numerics (no sampling)
same_inds_prob = zeros(length(k_R_vec), length(N_vec));

%h_x = h_x_vec(1); % TEMP WRONG !!!
for i=1:length(k_R_vec) % loop on different relatives relations
    k_R = k_R_vec(i);  % compute joint distribution of two family members: siblings/DZ and MZ twins    
    switch compute_mode
        case {'simulation', 'simulations'}
            qtl_vec = zeros(length(N_vec),iters); qtl_vec_R = zeros(length(N_vec),iters);
            inds_vec = zeros(length(N_vec),iters); inds_vec_R = zeros(length(N_vec),iters);
            for j=1:length(N_vec)
                data = zeros(iters,N_vec(j)); data_R = zeros(iters,N_vec(j)); % simulate family members

            for k=1:N_vec(j)
                tmp_data = simulate_correlated_gaussians(iters, 2, k_R*h_x_vec(j) + h_shared_env(j)); % correlation of sibs 
                data(:,k) = tmp_data(:,1); data_R(:,k) = tmp_data(:,2);
            end
            [qtl_vec(j,:) inds_vec(j,:)] = ...
                qtl_operation_internal(data, operation, operation_param, N_vec(j));
            [qtl_vec_R(j,:) inds_vec_R(j,:)] = ...
                qtl_operation_internal(data_R, operation, operation_param, N_vec(j));
            end
            same_inds_prob(i,:) = mean(inds_vec == inds_vec_R,2);
    end
    for j=1:length(N_vec)
        switch compute_mode
            case {'simulation', 'simulations'}
                qtl_R(i,j) = corr(vec2column(qtl_vec(j,:)), vec2column(qtl_vec_R(j,:))); % correlation between family members. That's not really lambda!
                
            case {'semi-analytic', 'numeric'}  % New: compute semi-analytic integration. Works only for max/min
                rho = k_R * h_x_vec(j) + 1*h_shared_env(j); % correlation includes shared environment 
                [mu_z sigma_z] = maxnormstat(N_vec(j));

                if(rho>1)
                    error('rho>1 !!!'), return;
                end
                if(rho==1) % degnerate case
                    qtl_R(i,j) = 1;
                else
                    y_max = mu_z + 5.*sigma_z; y_min = mu_z - 5.*sigma_z; % take 5 st.ds. above mean for MAX-of-Gaussians
                    qtl_R(i,j) = quad2d(@(z,z_R) z.*z_R.* ...
                        joint_two_gaussian_maxima_integrand_internal(z, z_R, N_vec(j), rho), ...
                        y_min, y_max, y_min, y_max);
                    qtl_R(i,j) = (qtl_R(i,j) - mu_z^2) / sigma_z^2; % normalize to get correlation
                end
        end
        
        %        should_be_one_normalized = quad2d(@(z,z_R) joint_two_gaussian_maxima_integrand_internal(z, z_R, N_vec(j), rho), ...
        %            -10, 10, -10, 10);
        
        
        %         x_grid = -5:0.1:5;
        %         x_mat = repmat(x_grid, length(x_grid), 1);
        %         y_mat = joint_two_gaussian_maxima_integrand_internal(x_mat, x_mat', N_vec(j), rho);
        %         figure; surf(x_mat, x_mat', y_mat);
        
        %        should_be_zero_symmetric = ...
        %            joint_two_gaussian_maxima_integrand_internal(2, 1, N_vec(j), rho) - ...
        %            joint_two_gaussian_maxima_integrand_internal(1, 2, N_vec(j), rho)
        %       should_be_zero_symmetric_cumulative = ...
        %            joint_two_gaussian_maxima_cumulative_internal(2, 1, N_vec(j), rho) - ...
        %            joint_two_gaussian_maxima_cumulative_internal(1, 2, N_vec(j), rho)
        
        
        %        if(abs(should_be_zero_symmetric) > 0.00001)
        %            xxx = 1241234
        %        end
        % %         %         joint_two_gaussian_maxima_integrand_internal([1 4], [2 5], N_vec(i), rho)
        % %         %         joint_two_gaussian_maxima_integrand_internal([1 4; 3 2], [2 5; 1 7], N_vec(i), rho)
        
    end
    if(i == -1) % plot
        figure; hold on;
        plot(qtl_vec, qtl_vec_R, '.');
        title(['QTL Corr of k_R=' num2str(k_R) ' is: ' num2str(qtl_R(i))]);
        xlabel('QTL individual'); ylabel('QTL family member');
    end
end

h_x_is = h_x_vec;
%%h_pop = 2 .* (qtl_R(1,:) - qtl_R(2,:)); % compute heritability using simple Falconer's estimator
switch epi_str
    case 'MZ'
        h_pop = qtl_R(1,:); % use just MZ twins !!!!
    case 'ACE'
        h_pop = 2*(qtl_R(1,:)-qtl_R(2,:)); % subtract DZ from MZ (account for common environment)       
    case 'ADE'
        h_pop = 4*qtl_R(2,:)-qtl_R(1,:); % subtract MZ from DZ (account for dominance)       
    case 'DZ'
        h_pop = qtl_R(2,:); % report just DZ twins !!!! that's not really a meaningfull estimate
    case 'PO' % correlate offpring with mid-parent 
        h_pop = 2*qtl_R(2,:); % assume no dominance (so parent-offspring is like sibling). What about shared envirounment???  
end
sprintf('N=%ld,h_x(one liab)=%0.2f%%, h_pop(trait)=%0.2f%%\n', N, 100*h_x_vec(1), 100*h_pop(1)); 
%h_pop = qtl_R(1,:); % Just use MZ twins ...



% Internal function: the joint density of two maximas of N Gaussians
% with the Gaussians correlated in pairs
%
% Input:
% x - value of first maximum
% y - value of second maximum
% N - number of Gaussians in each maximum
% rho - correlation coefficient between pairs of Gaussians
%
function ret = joint_two_gaussian_maxima_integrand_internal(z, z_R, N, rho)

m = size(z_R, 1); n = size(z_R, 2);
ttt = cputime;
if(rho > 1)
    ret = 0; error('rho>1'), return;
end
if(rho == 1) % degenerate case: perfect correlation 
    ret = zeros(size(z)); ret(z == z_R) = 9999999999999; return;
end
temp_int = reshape(mvncdf([z(:) z_R(:)], [0 0], [1 rho; rho 1]), m, n);
joint_max_time = cputime - ttt;
ret2 = repmat(temp_int.^(N-2), 1, 1) ./ (normcdf(z).^(N-1) .* sqrt(1-rho.^2)); % repmat (size(z_R, 1), 1)
z_R_normalized = (z_R-z.*rho)./sqrt(1-rho.^2); % normalize z_R conditioned on z
ret2 = ret2 .* (normpdf(z_R_normalized) .* repmat(temp_int, 1, 1) + ...
    (N-1) .* normcdf(z_R_normalized) .* ...
    exp(-z_R.^2./2) .* sqrt((1-rho.^2)./(2*pi)) .* normcdf((z-z_R.*rho)./sqrt(1-rho.^2)));
ret = maxnormpdf(z,N) .* ret2; %  .* z .* z_R; % multiplication by z*z_R should keep symmetry
% if(any(isnan(ret)))
%     wtf_isnan = 111
% end

% ret = N.*normpdf(x) .* normcdf(x).*(N-1); % compute g_N(x)
% Missing is the computation of g_N(x|y)
%ret2 = quadl(

% % %  %New approach: build BIG matrix!!!
% % % C = 0.5.*eye(2*N);
% % % for i=1:N
% % %     C(2*i-1,2*i) = rho;
% % % end
% % % C = C+C';
% % % for i=1:m
% % %     for j=1:n
% % %         ret(i,j) =



% temp_int = zeros(m,n); temp_int2 = zeros(m,n);
% ttt = cputime;
% for i=1:m % compute integral one by one
%     %    integrate_i = i
%     for j=1:n
% %        temp_int(i,j) = quadl(@(x)  norm_density_times_cumulative_internal(x, z_R(i,j), rho) , -10, z(i,j));        % That's the slowest part ...
%         temp_int(i,j) = mvncdf([z(i,j) z_R(i,j)]', [0 0]', [1 rho; rho 1]);
%     end
% end
% joint_max_time = cputime - ttt



% % Try again with cumsum cumsum2d(a)
% ttt2 = cputime;
% for i=1:m % compute integral one by one
% %        integrate_i2 = i
%     for j=1:n
%         temp_int2(i,j) = mvnpdf([z(i,j) z_R(i,j)]', [0 0]', [1 rho; rho 1]);
%     end
% end
% temp_int2 = cumsum2d(temp_int2);
% joint_numeric = cputime - ttt2
%
% max_diff = max(abs(temp_int(:) - temp_int2(:)))


% Internal integrad function
function ret = norm_density_times_cumulative_internal(x, z_R, rho)

ret = normpdf(x) .* normcdf((z_R - x.*rho)./sqrt(1-rho.^2));



