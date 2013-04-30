% Here we want to find the asymptotics of f in several parameters
% For now we do it for alpha -> 0 and alpha -> 1

% first take alpha to zero
TOL = 0.00000000001;
res = 0.00001; sigma = 1;
large_alpha_flag = 0; num_vars =500;
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2;

% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0; 
% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0; 

sigma_vec = [10.1:0.1:10.2];

a_vec = zeros(1,length(sigma_vec));
b_vec = zeros(1,length(sigma_vec));

alpha_vec = [res:res:20*res];
if(large_alpha_flag)
    alpha_vec = 1-alpha_vec;
end

figure; hold on;
iter=1;
for sigma=sigma_vec

    
    c_alpha_vec = norminv(1-0.5.*alpha_vec);
    x_alpha_vec = sqrt(1+sigma^2) .* c_alpha_vec;

    c_alpha_asym_vec = sqrt(-2 .* log(alpha_vec) - log(-2 .* log(alpha_vec)) +1* log(2/pi));
    % Give the asymptotic approximation for alpha -> 0 from Liat
    sig_minus = 1 - sqrt(1+sigma^2);
    sig_plus = 1 + sqrt(1+sigma^2);
    f_asym_vec = (-sigma / sig_minus) .* alpha_vec .^ (1 + 2*sig_minus/(sigma^2)) .* ...
        (sqrt(2/pi)./c_alpha_vec) .^ (-2*sig_minus/(sigma^2));

    % Now add the 'plus' part 
    f_asym_vec =  f_asym_vec + 0*(sigma / sig_plus) .* alpha_vec .^ (1 + 2*sig_plus/(sigma^2)) .* ...
        (sqrt(2/pi)./c_alpha_vec) .^ (-2*sig_plus/(sigma^2));
    
  %%  f_asym_vec = f_asym_vec ./ (2 + sig_minus.^2 ./ sigma.^2); % .* -(sigma^2) ./ (2 .* (sigma.^2 - sqrt(1+sigma.^2)));
          %%%f_asym_zuk_vec =  (sqrt(2/pi) * sigma ./ ( -sig_minus  )) .*
          %%%...
%    f_asym_zuk_vec = (sqrt(2/pi) * sigma^3 ./ ( -sig_minus .*  (sigma.^2+sig_minus.^2) )) .* ... 
%        alpha_vec .^ ((sig_minus./sigma).^2) .* (pi .* log(1./alpha_vec)).^ (sig_minus.^2./2.*sigma.^2) ./ c_alpha_vec;
        
   f_asym_zuk_vec = sqrt(2/pi) .* (-sigma ./ (sig_minus .* c_alpha_vec)) .* exp( -c_alpha_vec.^2 .* sig_minus .^ 2 ./ (2 .* sigma^2)   ) + ...
       sqrt(2/pi) .* (sigma ./ (sig_plus .* c_alpha_vec)) .* exp( -c_alpha_vec.^2 .* sig_plus .^ 2 ./ (2 .* sigma^2)   );
    
    f_asym_zuk_vec2 = 4 - 2*normcdf(-c_alpha_vec*sig_minus/sigma) - 2*normcdf(c_alpha_vec*sig_plus/sigma);
    
    
    
 %   f_asym_zuk_vec = (2-2*normcdf(-c_alpha_vec*sig_minus)); %%% ./ (1 + sig_minus.^2 ./ 1); % + ...
%     f_asym_vec = (2-2*normcdf(-c_alpha_vec*sig_minus)) .* (sigma^2) ./ (2 .* (1+sigma.^2 - sqrt(1+sigma.^2))) + ...
%                  (2-2*normcdf(c_alpha_vec*sig_plus)) .* (sigma^2) ./ (2 .* (1+sigma.^2 + sqrt(1+sigma.^2))) ;

    
    
%        (2-2*normcdf(c_alpha_vec*sig_plus)) .* (sigma^2) ./ (2 .* (1+sigma.^2 + sqrt(1+sigma.^2)));
    
    % Correction due to differentiation
%%%    f_asym_zuk_vec = f_asym_zuk_vec .* sigma^2 ./ (2 .* (sigma.^2 + sig_minus));
    
    
%     f_asym_zuk_vec_2 = sqrt(2/pi) .* (sigma ./ (sig_plus .* c_alpha_vec)) .* exp( -c_alpha_vec.^2 .* sig_plus .^ 2 ./ (2 .* sigma^2)   ); 
%      f_asym_zuk_vec = f_asym_zuk_vec_2 .* sigma^2 ./ (sigma.^2 + sig_plus.^2);
%     
%      f_asym_zuk_vec = f_asym_zuk_vec + f_asym_zuk_vec_2;
    
    f_vec = zeros(1, length(alpha_vec)); f_vec_sampled = zeros(1, length(alpha_vec));

    
    nsamples = 1/sigma^2;
    for i = 1:length(alpha_vec)
        int_right_limit = max(10, x_alpha_vec(i)*10);
        f_vec(i) = 1 - 2*quadl('f_int', c_alpha_vec(i), int_right_limit, TOL, [], x_alpha_vec(i), sigma)/alpha_vec(i);
        
        
          [kept_frac_dist f_all_genes f_all_genes_one_NTOP]= sample_kept_fraction_distribution(GAUSSIAN, TWO_SIDES,  TRUE_AND_SAMPLED, 1, nsamples, num_vars, alpha_vec(i), 20, R, FALSE);
          f_vec_sampled(i) = mean(kept_frac_dist)/(alpha_vec(i)*num_vars);
    end

    % Take the difference from one
    if(large_alpha_flag)
        f_vec = 1-f_vec;
    end

    
    
        rat_vec = f_asym_zuk_vec ./ f_vec;

    
    % Now plot the results :
    hold on;
    if(iter == 1)
%        plot((1-alpha_vec), (f_vec),'+');
         plot((alpha_vec), (f_vec),'+');
         plot(alpha_vec, f_asym_zuk_vec, 'r');
         plot(alpha_vec, f_asym_zuk_vec2, 'm*');
         plot(alpha_vec, alpha_vec, 'k');
         plot(alpha_vec, f_vec_sampled, 'g*');
         legend('numeric', 'asymptotic', 'new asym', 'alpha', 'sampled');
             xlabel('alpha'); ylabel('mean f');
    title(['loglog f in the alpha->0 regime for gaussian q for sigma = ' num2str(sigma)]);

%    else
%        plot((1-alpha_vec), (f_vec),'r.');
    end


    
    
    f_diff_vec = f_vec(2:end) + alpha_vec(2:end) .* (f_vec(2:end)-f_vec(1:end-1)) ./ res; 
    
    if(iter == 1)
        % Do another plot, here of f + alpha f'
        figure; hold on; 
        plot(alpha_vec(2:end), f_diff_vec, '+');
        plot(alpha_vec, f_asym_zuk_vec2, 'r');
        legend('numeric', 'asymptotic');
        xlabel('alpha'); ylabel('f + \alpha df/ d\alpha');
        title(['f + \alpha df/ d\alpha in the \alpha->0 regime for gaussian q for sigma = ' num2str(sigma)]);
    end
    
    
    
    %%%%%coef = fit(log(1-alpha_vec'), log(f_vec'),'a*x+b');
   
    
    %coef.a
    %coef.b
    %plot(log(1-alpha_vec'), coef.a.*log(1-alpha_vec')+coef.b, 'r');
    %figure; loglog(alpha_vec, f_vec); xlabel('alpha'); ylabel('mean f'); title('f in the alpha->0 regime for gaussian q loglog plot');

%%%%%    a_vec(iter) = coef.a;
%%%%%    b_vec(iter) = coef.b;
    iter=iter+1


end

return;


% Now plot a and b as function of sigma
figure; hold on; plot(sigma_vec, a_vec);plot(sigma_vec, exp(b_vec), 'r');
legend('a', 'b'); xlabel('sigma'); title('linear coefficients for f'); ylabel('coefs');

figure; hold on; plot(sigma_vec, exp(b_vec), 'r');
xlabel('sigma'); title('linear coefficient d for f'); ylabel('coef');