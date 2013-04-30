% Plot the predicted analyticity region of HMP entropy 

res = 0.001;
eps_vec = [0:res:1];
p_vec = [0:res:1];

zero_vec = zeros(1, length(eps_vec));
one_vec = ones(1, length(eps_vec));

figure; subplot(2,1,1); hold on; 
subplot(2,1,2); hold on; 
 plot(p_vec, zero_vec, 'r'); plot(p_vec, one_vec, 'r'); line([0.5,0.5],[0,1]); xlabel('p'); ylabel('\epsilon'); colorbar;
 
 
 % Now do the radius of convergence plot:
 rho_vec = 0.5 - (1-p_vec) ./ (2.*sqrt(1-2.*p_vec)); 
 
 subplot(2,1,1); hold on; plot(p_vec(1:495), -rho_vec(1:495), 'k'); plot(p_vec(end:-1:end-494), -rho_vec(1:495), 'k');
 xlabel('p'); ylabel('Estimated Convergence Radius');
 title('Estimated Convergence Radius');
 

