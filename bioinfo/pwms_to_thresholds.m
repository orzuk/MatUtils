% Compute a 'hit' threshold for each pwm
% such that on average,
% a fraction of alpha of random sequences (currently i.i.d.)
% will recieve a pwm score abouve this threshold
%
% Input:
% pwms - the pwms
% alpha - the fraction of sites passing the threshold.
% This flag is optional. If it isn't given, this number is computed from the pwms,
% based on its information content, via the formula: alpha = 1 / 2^IC
% iters - how many iterations (random sequences generated) to use
% background_model - a singletons model for generating background sequences (optional - default is [0.25 0.25 0.25 0.25])
%
% Output:
% pwms_thresholds - a threshold for each pwm such that a fraction alpha is above it
% pwms_alpha_vec - the fraction of sites expected to be above the thresholds
%
function [pwms_thresholds pwms_alpha_vec] = pwms_to_thresholds(pwms, alpha, iters, background_model, varargin)

if(~exist('iters', 'var') || isempty(iters))
    iters = 50000;
end
num_pwms = length(pwms);

if(~exist('alpha', 'var'))
    alpha = [];
end
if(~isempty(alpha))
    if(length(alpha) == 1)
        pwms_alpha_vec = repmat(alpha, num_pwms, 1);
    else
        pwms_alpha_vec = alpha;
    end
else
    pwms_alpha_vec = 1 ./ 2.^pwms_information_content(pwms);
end

L = 0;
for i=1:num_pwms
    L = max(L, length(pwms{i}));
end

if(~exist('background_model', 'var'))
    rand_seqs = ceil(rand(iters, L) .* 4);
else
    rand_seqs = weighted_rand(background_model, iters*L); rand_seqs = reshape(rand_seqs, iters, L);
    %    bp_stats = basecount(rand_seqs)
end

packed_rand_seqs = pack_seqs(rand_seqs);

pwms_thresholds = zeros(num_pwms, 1);
% % pwms_thresholds_gauss = zeros(num_pwms,1);  % temporary gaussian thresholds
at_home=0; % set the place to run
for i=1:num_pwms
    pwms{i} = single(pwms{i}); % force singles. doubles  do not yet work in the c function for some reason ...
    is_double = isa(pwms{i}, 'double');
    %%    mean_vec(i) = sum(sum(log(pwms{i}))) / 4;
    if(at_home) % didn't compile the 5-flags yet - this shouldn't be trusted
        loc_scores_mat = calc_all_sites_scores(double(log(pwms{i})), ...
            packed_rand_seqs, L, -10); % assume they're all equal and take the first
    else % at the broad. Seems to work, but only for singles ..(need to check again and see that scores match)
        loc_scores_mat = calc_all_sites_scores(log(pwms{i}), ...
            packed_rand_seqs, L, -10, is_double); % assume they're all equal and take the first
    end
    pwms_thresholds(i) = quantile(loc_scores_mat(:,1), 1-pwms_alpha_vec(i));

    % Another method: Just use Gaussian approximation. This should give a
    % (very) rough estimate of the threshold
    % %   pwm_score_mean = sum( mean  ( log(pwms{i}) ));
    % % pwm_score_var = sum( var  ( log(pwms{i}),1 ));
    % % pwms_thresholds_gauss(i) = norminv(1-pwms_alpha_vec(i)) * sqrt(pwm_score_var) + pwm_score_mean;
end

% t_gauss = pwms_thresholds_gauss



