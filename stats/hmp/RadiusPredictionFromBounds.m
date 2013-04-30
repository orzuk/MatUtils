% Predict entropy convergence radius from upper and lower bounds. 
% Note that this is still not rigrously proven to be working


tol = 0.000000000001; % Tolerance for being real

eps_start = -0.3; eps_end = 0.0-tol;
p_start = 0.391; p_end = 0.41-tol;

eps_vec = [eps_start:0.01:eps_end];
p_vec = [p_start:0.01:p_end];

C_N = zeros(length(p_vec), length(eps_vec));
cx_N = zeros(length(p_vec), length(eps_vec));

N=3;
i=1;
for eps = eps_vec
    C_N(:,i) = HMP_entropy_finite(eps, p_vec, N+1) - HMP_entropy_finite(eps, p_vec, N);
    cx_N(:,i) = HMP_entropy_finite_X1(eps, p_vec, N+1) - HMP_entropy_finite_X1(eps, p_vec, N);
    i = i+1;
end



C_N = real(C_N .* (C_N == abs(C_N)) );
cx_N = real(cx_N .* (abs(imag(cx_N)) < tol));


% Distinguish real zero from close to zero
C_N ( C_N == 0) = -1; 
cx_N ( cx_N == 0) = -1; 




% Try to calculate the radius analytically for N=infinity
%%%inf_rad_plus = 0.5 + sqrt(1+2.*p_vec.^2-3.*p_vec) ./ (2 .* (-1+2.* p_vec)); 
inf_rad_plus_eq_lam = 0.5 + sqrt(1-4.*p_vec+6.*p_vec.^2-4.*p_vec.^3) ./ (2 .* (-1+2.* p_vec)); 
inf_rad_minus = 0.5 - sqrt(1+2.*p_vec.^2-3.*p_vec) ./ (2 .* (-1+2.* p_vec)); 

inf_rad_plus = 0.5 - (1-p_vec)./(2.*sqrt(1-2.*p_vec));

figure; imagesc(C_N'); colorbar; ylabel(['epsilon = ' num2str(eps_start) ' - ' num2str(eps_end) ]); 
xlabel(['p = ' num2str(p_start) ' - ' num2str(p_end)]); title(['Upper N=' num2str(N)]);

figure; imagesc(cx_N'); colorbar; ylabel(['epsilon = ' num2str(eps_start) ' - ' num2str(eps_end) ]); 
xlabel(['p = ' num2str(p_start) ' - ' num2str(p_end)]); title(['Lower N=' num2str(N)]);

% Now try doing a decent plot 

C_N_pos = (C_N >= 0); C_N_pos(:,end) = 1;
cx_N_pos = (cx_N >= 0); cx_N_pos(:,end) = 1;

[C_N_val C_N_ind ] = max(C_N_pos');
[cx_N_val cx_N_ind ] = max(cx_N_pos');
figure; subplot(2,1,1); hold on; subplot(2,1,1); hold on; plot(p_vec, abs(eps_vec(C_N_ind)), 'b'); plot(p_vec, abs(eps_vec(cx_N_ind)), 'r'); 
plot(p_vec, min(666,p_vec ./ (1-2*p_vec)), 'g'); % plot(p_vec, HMMFittedRadius, 'm');
title(['Radius of Upperbound and LowerBound N=' num2str(N)]); legend('Upper', 'Lower', 'IID'); xlabel('p'); ylabel('radius');
subplot(2,1,2); hold on; plot(p_vec, abs(eps_vec(C_N_ind)), 'b'); plot(p_vec, abs(eps_vec(cx_N_ind)), 'r'); 
plot(p_vec,-inf_rad_plus, 'm'); plot(p_vec,-inf_rad_plus_eq_lam, 'c'); % plot(p_vec, -inf_rad_minus, 'k'); 
% plot(p_vec, HMMFittedRadius, 'm');
title(['Radius of Upperbound and LowerBound N=' num2str(N)]); legend('Upper', 'Lower', 'inf plus', 'inf plus eq lam'); xlabel('p'); ylabel('radius');