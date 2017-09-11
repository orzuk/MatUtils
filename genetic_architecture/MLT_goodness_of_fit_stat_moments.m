% Compute the difference in test statistic under the LT and MLT models.
% The inputs mu and h_x are ALWAYS the epidemiological disease parameters
% (not the single liability parameters) 
function [gof_mu_diff gof_sigma_diff pow] = ...
    MLT_goodness_of_fit_stat_moments(h_x, mu, N, K, n, alpha, true_model_str, ...
    iters, max_S)


AssignGeneralConstants();

% The true model index is always 1
h_x(2) = heritability_scale_change_MLT(N*h_x, K, N, mu, 'MLT'); % Get the heritabiity of each liability in MLT model
mu(2) = fminbnd(@(x) abs(binocdf(K-1, N, x)-(1-mu)), 0, 1); % find mu_l that keeps the prevalence
N = [1 N];
switch true_model_str
    case 'LT'
        LT_ind=1; model_str = {'LT', 'MLT'};
    case 'MLT'
        h_x = h_x(2:-1:1); mu = mu(2:-1:1); N = N(2:-1:1);
        LT_ind=2; model_str = {'MLT', 'LT'};
        
        %
        %         h_x(2) = heritability_scale_change_MLT(N*h_x, K, N, mu, 'LT'); % Get the heritabiity of each liability in MLT model
        %         mu(2) = 1-binocdf(K-1, N, mu);  % find mu_l that keeps the prevalence
end
MLT_ind = 3-LT_ind; 
x_mu(1) = norminv(1-mu(1));
x_mu(2) = norminv(1-mu(2)); % set threshold


% mu = int(phi(x) * (1-2p(x))*(q(x)-p(x))/ (p(x)*(1-p(x))
x_min = -5; x_max = 5; % integral limits (in standard deviations)

const_p_flag=0;
if(const_p_flag)% Do some tests:
    prevalence_LT = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
        z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu), ...
        x_min, x_max, x_min, x_max);
    prevalence_MLT = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
        z_expected_given_two_x_MLT(x1, x2, 2, h_x_l, mu_l, x_mu_l), ...
        x_min, x_max, x_min, x_max);
    
    
    gof_mu_LT = 1;
    gof_mu_MLT = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
        ((1-2*z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) .* ...
        z_expected_given_two_x_MLT(x1, x2, N, h_x_l, mu_l, x_mu_l) + ...
        z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu).^2) ./ ...
        ( z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu) .* ...
        (1-z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) ), ...
        x_min, x_max, x_min, x_max); % this gives you the DIFFERENCE in mu!!!
    
    gof_sigma_LT = 1 - 2*quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ( ...
        z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)), ...
        x_min, x_max, x_min, x_max); % This is zero!! (since we average over p)
    gof_sigma_MLT = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
        (1-2*z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) .* ...
        z_expected_given_two_x_MLT(x1, x2, N, h_x_l, mu_l, x_mu_l) .* ...
        (1-z_expected_given_two_x_MLT(x1, x2, N, h_x_l, mu_l, x_mu_l)) ./ ...
        ( z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu) .* ...
        (1-z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) ), ...
        x_min, x_max, x_min, x_max);
    
    
    gof_mu_diff = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
        (1-2*z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) .* ...
        (z_expected_given_two_x_MLT(x1, x2, N, h_x_l, mu_l, x_mu_l) - ...
        z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) ./ ...
        ( z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu) .* ...
        (1-z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) ), ...
        x_min, x_max, x_min, x_max); % this gives you the DIFFERENCE in mu!!!
    gof_sigma_diff = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ( ...
        2.*z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)-1 + ...
        (1-2*z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) .* ...
        z_expected_given_two_x_MLT(x1, x2, N, h_x_l, mu_l, x_mu_l) .* ...
        (1- z_expected_given_two_x_MLT(x1, x2, N, h_x_l, mu_l, x_mu_l)) ./ ...
        ( z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu) .* ...
        (1-z_expected_given_two_x_MLT(x1, x2, 1, h_x, mu, x_mu)) ) ), ...
        x_min, x_max, x_min, x_max); % this gives you the DIFFERENCE in sigma!!!
    
else % New computation: integrate also over p and q - we need also the test
    % statistic to include that (?)
    % mu + phi(x)*q(x)*(1-2p(x))/p(x)
    
    gof_mu_numeric = zeros(2); gof_sigma_numeric = zeros(2);
    for i=1:2 % loop on two different models GENERATING data (LT and MLT)
        for j=1:2  % loop on two different models ASSUMED AS NULL data (LT and MLT)
            gof_mu_numeric(i,j) = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
                (1-2*z_expected_given_two_x_MLT(x1, x2, N(j), h_x(j), mu(j), x_mu(j))) .* ...
                z_expected_given_two_x_MLT(x1, x2, N(i), h_x(i), mu(i), x_mu(i)) ./ ...
                z_expected_given_two_x_MLT(x1, x2, N(j), h_x(j), mu(j), x_mu(j)), ...
                x_min, x_max, x_min, x_max) + mu(LT_ind);
            gof_sigma_numeric(i,j) = sqrt( quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
                ( z_expected_given_two_x_MLT(x1, x2, N(j), h_x(j), mu(j), x_mu(j)).^2 - ...
                4.*z_expected_given_two_x_MLT(x1, x2, N(j), h_x(j), mu(j), x_mu(j)).* ...
                z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)) + ...
                (1-2*z_expected_given_two_x_MLT(x1, x2, N(j), h_x(j), mu(j), x_mu(j))).^2 .* ...
                z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)) ./ ...
                z_expected_given_two_x_MLT(x1, x2, N(j), h_x(j), mu(j), x_mu(j)).^2 ), ...
                x_min, x_max, x_min, x_max) +2*mu(LT_ind) - gof_mu_numeric(i,j)^2 );
            
            
            %     gof_mu_MLT = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
            %         (1-2*z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1))) .* ...
            %         z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)) ./ ...
            %         z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)), ...
            %         x_min, x_max, x_min, x_max) + mu(LT_ind);
            %     gof_mu_diff = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
            %         (1-2*z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1))) .* ...
            %         z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)) ./ ...
            %         z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)), ...
            %         x_min, x_max, x_min, x_max) + 2*mu(LT_ind)-1; % This is MLT - LT
            
            %     gof_sigma_LT = sqrt( quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
            %         ( -3*z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)).^2 + ...
            %         (1-2*z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1))).^2 ./ ...
            %         z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)) ), ...
            %         x_min, x_max, x_min, x_max) -1-mu(LT_ind)^2+4*mu(LT_ind) );
            %
            %
            %     gof_sigma_diff = quad2d(@(x1,x2) normpdf(x1) .* normpdf(x2) .* ...
            %         ( 4.*z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)).* ...
            %         (z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)) - ...
            %         z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2))) + ...
            %         (1-2*z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1))).^2 .* ...
            %         (z_expected_given_two_x_MLT(x1, x2, N(2), h_x(2), mu(2), x_mu(2)) - ...
            %         z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1))) ./ ...
            %         z_expected_given_two_x_MLT(x1, x2, N(1), h_x(1), mu(1), x_mu(1)).^2 ), ...
            %         x_min, x_max, x_min, x_max);
        end % loop on data generating models
    end % loop on data analyzing model
    
    rand_str = 'grid'; 
    gof_mu_diff = gof_mu_numeric(MLT_ind,LT_ind) - gof_mu_numeric(LT_ind,LT_ind);
    gof_sigma_diff = gof_sigma_numeric(1,2) - gof_sigma_numeric(1,1);
%    iters = 100^2; % 1000000; % Compute moments also empirically
    switch rand_str
        case 'rand'
    X = randn(iters,2); 
        case 'grid'
            [xx yy] = meshgrid(norminv(((1:sqrt(iters))-0.5) ./ sqrt(iters)));
            X = [mat2vec(xx) mat2vec(yy)];     
    end
    z_expected = zeros(iters,2); z_observed = zeros(iters,2);
    gof_mu_empirical = zeros(2); gof_sigma_empirical = zeros(2);
    
    for i=1:2 % loop on two different models GENERATING data (LT and MLT)
        if(N(i) >= 2) % loci are in different liabilities
            z_expected(:,i) =  1 - ((1-mu(i)).^(N(i)-2) .* ...
                normcdf( (x_mu(i) - X(:,1) .* sqrt(h_x(i))) ./ sqrt(1-h_x(i))) .* ...
                normcdf( (x_mu(i) - X(:,2) .* sqrt(h_x(i))) ./ sqrt(1-h_x(i))));  % z_expected
        else
            z_expected(:,i) = 1 - (normcdf( (x_mu(i) - X(:,1) .* sqrt(h_x(i)) - ...
                X(:,2) .* sqrt(h_x(i))) ./ ...
                sqrt(1-2*h_x(i)))); % z_expected
        end
        z_observed(:,i) = rand(iters, 1) < z_expected(:,i); % simulate according to correct model
    end
    for i=1:2
        for j=1:2 % loop on two different models ASSUMED AS NULL data (LT and MLT)
            stat_vec{i,j} = (z_observed(:,i) - z_expected(:,j)).^2 ./ z_expected(:,j);
            gof_mu_empirical(i,j) = mean(stat_vec{i,j});
            gof_sigma_empirical(i,j) = std(stat_vec{i,j});
            stat_vec_replicas{i,j} = cumsum(stat_vec{i,j});
            stat_vec_replicas{i,j} = diff(stat_vec_replicas{i,j}(1000:1000:end));

        end
    end
    
    
    % Compare empirical to theoretical moments of test statistic
    for plot_ind = 1:2
        figure; hold on;
        
        if(plot_ind == 1)
            plot(gof_mu_numeric(:), gof_mu_empirical(:), '*');
            x_vec = linspace(min(min(gof_mu_numeric(:)),min(gof_mu_empirical(:))), ...
                max(max(gof_mu_numeric(:)),max(gof_mu_empirical(:))), 10);
            
        else
            plot(gof_sigma_numeric(:), gof_sigma_empirical(:), '*');
            x_vec = linspace(min(min(gof_sigma_numeric(:)),min(gof_sigma_empirical(:))), ...
                max(max(gof_sigma_numeric(:)),max(gof_sigma_empirical(:))), 10);
        end
        plot(x_vec, x_vec, 'k');
        
        for i=1:2
            for j=1:2
                if(plot_ind==1)
                    text(gof_mu_numeric(i,j)+0.001, gof_mu_empirical(i,j)+0.001, ...
                        ['(true:' model_str{i} ', null:' model_str{j} ')']);
                else
                    text(gof_sigma_numeric(i,j)+0.005, gof_sigma_empirical(i,j)+0.005, ...
                        ['(true:' model_str{i} ', null:' model_str{j} ')']);
                end
                
            end
        end
        
        %         plot(gof_mu_LT, gof_mu_LT_empirical, '*');
        %         plot(gof_mu_MLT, gof_mu_MLT_empirical, '*r');
        %         plot(gof_sigma_LT, gof_sigma_LT_empirical, '*m');
        %         plot(gof_sigma_MLT, gof_sigma_MLT_empirical, '*g');
        %        legend('mu-LT', 'mu-MLT', 'sigma-LT', 'sigma-MLT');
        xlabel('numeric'); ylabel('simulated');
        if(plot_ind == 1)
            title('\mu for LT vs. MLT empirical and theoretical');
        else
            title('\sigma for LT vs. MLT empirical and theoretical');
        end
        
        
    end % loop on plot ind
    figure; hold on; legend_vec = cell(1,4); % Plot empirical histogram of test statistic
    for i=1:2
        for j=1:2
                hist_density(stat_vec{i,j}, 500, color_vec(i*2+j-1));
                legend_vec{(i-1)*2+j} = ['(true:' model_str{i} ', null:' model_str{j} ')'];
                
                
        end
    end    
%    xlim([0 2]);
    legend(legend_vec); xlabel('S_{gof}'); ylabel('Freq.'); 
    title('Distribution of GOF test statistic under the different models'); 

    figure; hold on; legend_vec = cell(1,4); % Plot empirical histogram of test statistic
    for i=1:2
        for j=1:2
                hist_density(stat_vec_replicas{i,j}, 50, color_vec(i*2+j-1));
                legend_vec{(i-1)*2+j} = ['(true:' model_str{i} ', null:' model_str{j} ')'];
                
                
        end
    end    
%    xlim([0 20]);
    legend(legend_vec); xlabel('S_{gof}'); ylabel('Freq.'); 
    title('Distribution of GOF test statistic n=1000 under the different models'); 
    
    
end % if const_p_flag
    


% Compute power
pow=1;
effect_size_mu_diff = gof_mu_diff
pow = 1 - normcdf(  (gof_sigma_numeric(LT_ind,LT_ind).*norminv(1-alpha)-sqrt(n).*gof_mu_diff) ./ ...
    gof_sigma_numeric(MLT_ind,MLT_ind))

%
% N = 10; p=0.1;
% A = eye(N);
% for i=2:N-1
%     A(i,i-1) = p-1; A(i,i+1)=-p;
% end
% b = zeros(N,1); b(end) = 1;
% x = linsolve(A,b) % b' / A;
% A*x'-b
%
