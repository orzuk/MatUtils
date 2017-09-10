% Compute the joint distribution of two genotypes from relative risk 
% and risk-allele-frequency
%
% Input:
% f_vec - risk allele frequencies Pr(x=1)
% lambda_R_vec - familial relative risk Pr(z=1|x=1) / Pr(z=1|x=0)
% mu_vec - prevalence Pr(z=1)
%
% Output:
% p_x_x_R - joint proability of disease for two relative (row vector of size 4 for each input)
%
function p_x_x_R = lambda_r_to_p_z_x_marginal( lambda_R_vec, mu_vec)

p_x_x_R = zeros(length(mu_vec), 4);

mu_vec = vec2column(mu_vec);
lambda_R_vec = vec2column(lambda_R_vec); 

p_x_x_R(:,1) = 1-2.*mu_vec + mu_vec.^2 .* lambda_R_vec; % Pr(x=0,x_R=0)
p_x_x_R(:,2) = mu_vec .* (1-mu_vec.*lambda_R_vec); % Pr(x=0,x_R=1)
p_x_x_R(:,3) = mu_vec .* (1-mu_vec.*lambda_R_vec); % Pr(x=1,x_R=0)
p_x_x_R(:,4) = mu_vec.^2 .* lambda_R_vec; % Pr(x=1,x_R=1)

