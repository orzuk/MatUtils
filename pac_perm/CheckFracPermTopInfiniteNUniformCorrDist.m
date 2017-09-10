% A script for computing overlap distribution when correlations are uniformly distributed

% Here do only maple integral to show a pic.
res = 1; max_sigma = 10; alpha_res = 0.1;
alpha_vec = [alpha_res:alpha_res:0.5];
sigma_vec = [res:res:25]; %%% [res:res:max_sigma];



maple_int_res = zeros(length(alpha_vec),length(sigma_vec));

Xalpha = zeros(length(alpha_vec),length(sigma_vec)); 
for j=1:length(alpha_vec)
    j
    
    for i = 1:length(sigma_vec)
        
        % First determine Xalpha. It must be more than zero !!!! 
        Xalpha(j,i) = fsolve('JointDensAlphaFracUniformDist', 0,  optimset('fsolve'), sigma_vec(i), alpha_vec(j));
        JointDensAlphaFracUniformDist(Xalpha(j,i), sigma_vec(i), alpha_vec(j))
        
        maple_int_res(j,i) =  (1.0/alpha_vec(j))*JointDensFracUniformDist(Xalpha(j), sigma_vec(i), alpha_vec(j));
%%%%        maple_int_res(j,i) = (1.0/alpha_vec(j))*quad('JointDensFrac', C_alpha(j), 99999, [], [], sigma_vec(i), C_alpha(j));
    end
end

figure; hold on; %%%imagesc(  alpha_vec, sigma_vec,  maple_int_res); colorbar;  AXIS([alpha_res 0.5 res sigma_vec(end)]);
imagesc( sigma_vec,  alpha_vec, maple_int_res); colorbar;   AXIS([res sigma_vec(end) alpha_res 0.5 ]);
xlabel('Sigma'); ylabel('frac. alpha');  
title('Fraction kept Results for uniform correlations'); 


figure; hold on; 
imagesc( sigma_vec,  alpha_vec, Xalpha); colorbar;   AXIS([res sigma_vec(end) alpha_res 0.5 ]);
xlabel('Sigma'); ylabel('frac. alpha');  
title('Xalpha'); 
%%%%%%figure; hold on; plot(alpha_vec, Xalpha, '*'); title('Xalpha as a function of alpha'); 