% function test_goodness_of_fit()

p = [1/3 2/3]; % 'True' p we assume from which we sample data
q = [1/4 3/4]; % alternative q we actually draw data from
sigma = 0.0001; % st.d. of assumed vs. computed p
n = 100000;
iters = 1000;
stat_str = 'variance';

for i=1:iters
    if(mod(i,10)==0)
        run_iter =i
    end
    p_vec = p(1) + (rand(n,1) > 0.5) * (p(2)-p(1)); % expected
    q_vec = min(1, randn(n,1) * sigma + p_vec); % This is where we draw values from
    z_vec = rand(n,1)< q_vec;
    
    %    s_gof(i) = sum((z_vec - p_vec).^2 ./ p_vec);
    
    p_vec(:) = 0.5; % Change expected 
    switch stat_str
        case 'expected'
            s_gof(i) = sum((z_vec - p_vec).^2 ./ p_vec);
            stat_mu = n - sum(p_vec); % Compute statistic moment
            stat_sigma = sqrt(  sum((1-2.*p_vec).^2.*(1-p_vec) ./ p_vec) );
        case 'variance'
            s_gof(i) = sum((z_vec - p_vec).^2 ./ (p_vec.*(1-p_vec)));
            stat_mu = n; % Compute statistic moment
            stat_sigma = sqrt(n-  sum(p_vec) );
            
    end
    
    s_gof(i) =   (s_gof(i) - stat_mu) / stat_sigma;
end

figure; hist_density(s_gof, 50);
x_vec = -5:0.01:5; % :1:n_samples_vec(1);
plot(x_vec, normpdf(x_vec), 'r'); % , stat_mu, stat_sigma), 'r');
xlabel('goodness-of-fit-test-statistic'); ylabel('freq.');
title('Fit of wrong model');
legend('Observed', 'Guass. Approx.');   %            ['\chi^2_{' num2str(n_samples_vec(1)-1) '}']);

