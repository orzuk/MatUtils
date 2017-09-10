% Allele frequency distribution adjusted for power
%
% Input:
% x - allele frequency
% s - selection coefficient
% beta - effect size
% N - population size
% num_cases - number of individual cases
% num_controls - number of individual controls
% alpha - pval cutoff for declaring significance
% two_side_flag - take RAF or MAF
% scale_mode - linear or logarithmic scale 
% use_power - include power estimation in density
%
% Output:
% g - probability of allele frequency being x
%
function g = allele_freq_spectrum_power_corrected(x, s, N, beta, mu, ...
    num_cases, num_controls, alpha, two_side_flag, scale_mode, use_power)

if(~exist('two_side_flag', 'var') || isempty(two_side_flag))
    two_side_flag = 1; % default is MAF !!!
end
if(~exist('scale_mode', 'var') || isempty(scale_mode)) % set default as linear
    scale_mode = 'linear';
end
if(~exist('use_power', 'var') || isempty(use_power)) % set default as linear
    use_power = 1;
end
if(length(beta) == 1)
    beta = repmat(beta, length(x), 1);
end
if(length(s) == 1)
    s = repmat(s, length(x), 1);
end

trait_type = 'QTL'; num_points = length(x);
switch trait_type
    case 'QTL'
        p_mat = [vec2column(x) vec2column(beta) repmat([0 1], num_points, 1)]; % assume variance one
        test_type = 'single-locus';
        test_stat = 'chi-square-QTL-analytic';
    case 'binary'
        p_mat = genetic_relative_risk_to_p_z_x_marginal(x, beta, mu);
        test_type = 'armitage';
        test_stat = 'chi-square-analytic';
end

g = zeros(num_points,1); 
for i=1:num_points
    switch scale_mode
        case 'linear'
            g(i) = exp(allele_freq_spectrum(x(i), s(i), N, ...
                two_side_flag, 'log'));
            if(use_power)
                g(i) = g(i) .* compute_association_power(p_mat(i,:), ...
                    num_cases, num_controls, alpha, [], test_type, test_stat); % , sampling_type); % multiply spectrum by power
            end
        case 'log'
            g(i) = allele_freq_spectrum(x(i), s(i), N, ...
                two_side_flag, 'log');
            if(use_power)
                g(i) = g(i) +  ...
                    log(compute_association_power(p_mat(i,:), ...
                    num_cases, num_controls, alpha, [], test_type, test_stat)); % , sampling_type); % multiply spectrum by power
            end
    end % switch scale mode 
end % loop on points 


