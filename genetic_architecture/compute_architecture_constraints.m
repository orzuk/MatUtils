% Compute several quantities for a genetic architecture
% 
% Input:  
% coupling_str - type of coupling between distributions
% coupling_param - parameters for coupling 
% N - effective population size 
% s_vec - vector of possible selection coefficients (NEGATIVE for puryfying selection)
% beta_vec - vector of possible effect sizes 
% 
% Output: 
% mean_phenotype_change - change in phenotype  
% mean_square_phenotype_change - change in heritability 
% skewness - third moment
% kurtosis - fourth moment 
% intermediate_heritability - intermediate heritability 
%
function [mean_phenotype_change, mean_square_phenotype_change ...
    skewness, kurtosis, intermediate_heritability] = compute_architecture_constraints(coupling_str, coupling_param, N, s_vec, beta_vec)

s_min = min(s_vec); s_max = max(s_vec); beta_min = min(beta_vec); beta_max = max(beta_vec); S_vec = 4.*N.*s_vec;

[s_mat beta_mat] = meshgrid(s_vec, beta_vec);
Z = coupling_moment(s_mat, beta_mat, N, 'mean_phenotype');
% Y = g(s_mat, beta_mat, N, coupling_str, coupling_param); % Compute joint distribution of beta and s 
Y = joint_selection_effect_size_allele_freq(s_mat, beta_mat, ... % allele frequency distribution
    N, coupling_str, coupling_param); % effective population size
W = Y .* coupling_moment(s_mat, beta_mat, N, 'mean_phenotype');

figure; imagesc(s_vec, beta_vec, Y); colorbar; xlabel('s'); ylabel('\beta'); title('g(s,\beta)');
figure; imagesc(s_vec, beta_vec, Z); colorbar; xlabel('s'); ylabel('\beta'); title('other-stuff');
figure; imagesc(s_vec, beta_vec, W); colorbar; xlabel('s'); ylabel('\beta'); title('integrand');


mean_phenotype_change = quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s, beta, N, coupling_str, coupling_param) ...
    .* coupling_moment(s, beta, N, 'mean_phenotype'), s_min, s_max, beta_min, beta_max) ./ ...
    quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s,beta, N, coupling_str, coupling_param),...
    s_min, s_max, beta_min, beta_max); % normalize


mean_square_phenotype_change = quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s, beta, N, coupling_str, coupling_param) ...
    .* coupling_moment(s, beta, N, 'mean_square_phenotype'), s_min, s_max, beta_min, beta_max) ./ ...
    quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s,beta, N, coupling_str, coupling_param),...
    s_min, s_max, beta_min, beta_max); % normalize


intermediate_heritability = quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s, beta, N, coupling_str, coupling_param) ...
    .* coupling_moment(s, beta, N, 'intermediate_heritability'), s_min, s_max, beta_min, beta_max) ./ ...
    quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s,beta, N, coupling_str, coupling_param),...
    s_min, s_max, beta_min, beta_max); % normalize


skewness = quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s, beta, N, coupling_str, coupling_param) ...
    .* coupling_moment(s, beta, N, 'skewness'), s_min, s_max, beta_min, beta_max) ./ ...
    quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s,beta, N, coupling_str, coupling_param),...
    s_min, s_max, beta_min, beta_max); % normalize


kurtosis = quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s, beta, N, coupling_str, coupling_param) ...
    .* coupling_moment(s, beta, N, 'kurtosis'), s_min, s_max, beta_min, beta_max) ./ ...
    quad2d(@(s,beta) joint_selection_effect_size_allele_freq(s,beta, N, coupling_str, coupling_param),...
    s_min, s_max, beta_min, beta_max); % normalize


% % %         ;s
% % % g(
% % % ; s)
% % % 2e?S + S ? 2
% % % 1 ? e?S d
% % % ds = 0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Kernel function to multiply density g by to get constraint on architecture
% 
% Input: 
% s - selection coefficient
% beta - effect size
% N - effective population size 
% moment_str - type of constraint
%
% Output: 
% M - value of Kernel function 
% 
function M = coupling_moment(s, beta, N, moment_str)

if(isvector(beta))
    beta = repmat(beta, length(s), 1);
end

S = 4.*N.*s; % compute big S from little s 
[unique_s, I, J] = unique(s); % reduce number of computations
switch moment_str
    case 'mean_phenotype'  % K(s, \beta) = beta * (3e^(-S) + 2S - 3) / (1-e^(-S))
                           % New: should be: K(s, \beta) = beta * S / (1-e^(-S))
        integral_vec = absorption_time_by_selection(unique_s, 1, N, 0, 1, 'var'); % this gives int f(1-f) t_s(f) df 
        integral_vec = integral_vec(J);
        if(~isvector(s))
            integral_vec = reshape(integral_vec, size(s));
        end
        M = beta .* ( -2.*integral_vec - 1);
        %        (3.*exp(-S) + 2.*S - 3) ./ (1 - exp(-S));

    case 'mean_square_phenotype' % K(s, \beta) = beta^2 * (e^(-S) + S - 1) / (1-e^(-S))
        integral_vec = absorption_time_by_selection(unique_s, 1, N, 0, 1, 'var');
    
    
    case 'intermediate_heritability' % K(s, \beta) = beta^2 * (e^(-S) + S - 1) / (1-e^(-S)) - 1/2. Note: this is the same (up to a constant) as the previous one !!!! 
        integral_vec = absorption_time_by_selection(unique_s, 1, N, 0, 1, 'var');
        integral_vec = integral_vec(J);
        M = beta.^2 .* ( -2.*integral_vec - 1);
        
                
    case 'skewness'  % K(s, \beta) = beta^3 * ((S^2-1)(e^(-S) + S - 1)-S^2) / (S^2(1-e^(-S)))
        integral_vec = absorption_time_by_selection(unique_s, 1, N, 0, 1, 'var');
%        integral_vec = integral_vec(J);
    
    
    
    case 'kurtosis' % K(s, \beta) = beta^4 * ((S^3-3S+12)(e^(-S) + S - 1)-S^2(S+1)) / (S^3(1-e^(-S)))
        integral_vec = absorption_time_by_selection(unique_s, 1, N, 0, 1, 2);
        if(~isvector(s))
            integral_vec = reshape(integral_vec, size(s));
        end
        M = beta .* integral_vec;        
end

