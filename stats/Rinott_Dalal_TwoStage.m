% Analysis of constants in Rinot and Dalal procedures 

k_vec = 1:200; % list of k values 
p = 0.69; % probability of getting maximum correctly 
num_k = max(k_vec); 
nu = 2; % deg. freedom for t-distribution 

h_k_rinott = zeros(1, num_k+1);
h_k_dalal = zeros(1, num_k+1);

for k=1:length(k_vec)
    h_k_dalal(k+1) = fsolve(@(x) two_stage_integral_dalal(x, k, nu)-p, h_k_dalal(k)); 
    h_k_rinott(k+1) = fsolve(@(x) two_stage_integral_rinott(x, k, nu)-p, h_k_rinott(k));         
end
h_k_dalal = h_k_dalal(2:end);
h_k_rinott = h_k_rinott(2:end);


figure; plot(k_vec, h_k_dalal, '*'); 
hold on; plot(k_vec, h_k_rinott, 'r*'); 
legend({'dalal', 'rinott'}); legend('boxoff'); 
xlabel('k'); ylabel('h(k)'); 

figure; loglog(k_vec, h_k_dalal, '*'); 
hold on; loglog(k_vec, h_k_rinott, 'r*'); 
legend({'dalal', 'rinott'}); legend('boxoff'); 
xlabel('k'); ylabel('h(k)'); 

figure; plot(k_vec, h_k_rinott ./ h_k_dalal, 'g*'); 
xlabel('k'); ylabel('h_{rinnot}(k) / h_{dalal}(k)'); 

% Now implement procedures and call them: 

N0=100; PCS = 0.95; Delta = 0.1; n_pop = 5; mu_vec = zeros(1, n_pop); mu_vec(end) =  Delta; sigma_vec = ones(1, n_pop); iters = 1000; 
[max_I, PCS_D, h_D] = dalal_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'dalal');
[max_I_R, PCS_R, h_R] = dalal_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, 'rinott');
 
PCS_D
PCS_R

